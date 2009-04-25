#include <Driver.d/Sower.h>
#include <Utils.d/BinFileHandler.h>
#include <algorithm> // PJSA: for sgi intel

#ifdef SOWER_SURFS
#include <Utils.d/resize_array.h>
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#endif

//#define SOWER_DEBUG
//#define SOWER_DEBUG_SURFS

Sower::Sower(Connectivity* subToElem, Elemset& eset, int nClus, ResizeArray<SurfaceEntity*>* Surfs)
{
//Warning, subToElem is packed and eset is not ...
  subToClus = 0;
  tocRead=false;
#ifdef SOWER_DEBUG
  cerr << " ** Sower Constructor : " << endl;
  cerr << " ** subdomain to elements connectivity : " << endl;
  subToElem->print();
#endif
  Connectivity* eTos = subToElem->reverse();
#ifdef SOWER_DEBUG
  cerr << " ** elements to subdomains connectivity : " << endl;
  eTos->print();
#endif 
  Connectivity* sTos = subToElem->transcon(eTos);
#ifdef SOWER_DEBUG
  cerr << " ** subdomains to subdomains connectivity : " << endl;
  sTos->print();
#endif   
  Connectivity* eToN = new Connectivity(&eset);
#ifdef SOWER_DEBUG
  cerr << " ** elements to nodes connectivity : " << endl;
  eToN->print();
#endif   
  Connectivity* sToN = subToElem->transcon(eToN);
#ifdef SOWER_DEBUG
  cerr << " ** subdomains to nodes connectivity : " << endl;
  sToN->print();
#endif     
  Connectivity* nToS = sToN->reverse();
#ifdef SOWER_DEBUG
  cerr << " ** nodes to subdomain connectivity : " << endl;
  nToS->print();
#endif   
  // number of elements in a subdomain as weight
  int nSub = subToElem->csize();
  typedef long long superlong;
  long long * sizes = new superlong[nSub];
  for(int i = 0; i<nSub ; i++)
    sizes[i] = subToElem->num(i);
  
  clusToSub = greedy(sTos,nToS,sToN,sizes,nSub,nClus); 
  //addParentToChildData<Elemset*,Connectivity*>(ELEMENTS,CLUSTER,0,&eset,clusToSub);

#ifdef SOWER_DEBUG
  cerr << " ** cluster to subdomain connectivity : " << endl;
  clusToSub->print();
#endif   
  // PJSA: write clusToSub file
  BinFileHandler file("subdomains", "w");
  writeSubFile(file);
  // PJSA: write subToElem & subToNode file (temporary fix, replace with distributed decomposition)
  BinFileHandler file2("decomposition", "w");
  subToElem->write(file2);
  sToN->write(file2);
  // PJSA: write connectivity file
  BinFileHandler file3("connectivity", "w");
  clusToSub->write(file3);
  Connectivity *subToSub = sToN->transcon(nToS);
  subToSub->write(file3);
  delete subToSub;
  int nGlobNodes = nToS->csize();
  file3.write(&nGlobNodes,1);
  int nGlobElems = eTos->csize();
  file3.write(&nGlobElems,1);

  nCluster = clusToSub->csize(); /* true number of cluster */
  entries[ELEMENTS_TYPE] = new GenDataStruct<Elemset*,ElemsetIO>(&eset,/*clusToSub->transcon(*/subToElem/*)*/);

  numSubdomains = subToElem->csize();
  
  nSurfaces = 0;
  Surfaces  = Surfs;
  if(Surfaces) { // count real number of surfaces
    for(int i=0; i<Surfaces->max_size(); i++)
      if((*Surfaces)[i]) nSurfaces++;
  }
  if(nSurfaces) cerr << " ... "<<nSurfaces<<" surfaces have been defined in input file"<<endl;
  // Write (global) number of surfaces into connectivity files
  file3.write(&nSurfaces,1);

  delete eTos;
  delete sTos;
  delete eToN;
  delete sToN;
  delete nToS;
}

Sower::~Sower()
{
  for(map<TypeTag,DataStruct*>::iterator it = entries.begin(); it != entries.end(); it++)
    if((*it).second != 0)
      delete((*it).second);
  if(clusToSub)
    delete clusToSub;
}

void Sower::printDebug()
{
  std::cerr << " def : END=0,NODES=1,ELEMENTS=2,ATTRIBUTES=3,FORCES=4,MATERIALS=5,DISPLACEMENTS=6,\
	     EFRAMES=7,COMPOSITE=8,CFRAMES=9,CLUSTER=10 etc " << endl;
  for(map<TypeTag,DataStruct*>::iterator it = entries.begin(); it != entries.end(); it++)
    {
      std::cerr << endl << " Cluster to " << (*it).first << " connectivity : " << endl;
      (*it).second->getClusterToData(clusToSub)->print();
      std::cerr << endl << " SubDomain to " << (*it).first << " connectivity : " << endl;
      (*it).second->getSubToData()->print();
    }
}

void Sower::nomask(Connectivity* mn, int* maskno, int nelem) 
{
  int i;
  for(i = 0; i < mn->num(nelem); i++) maskno[(*mn)[nelem][i]]--;
}

Connectivity * 
Sower::greedy(Connectivity *etoe, Connectivity *ntoe,
              Connectivity *eton, long long *sizes, 
	      int decompNumSubdomains, int expectedNumSubdomains)  
{
  int i, minw, minw2;
  int node, node2, ntot, nelem, ptr1, ptr2;
  int *maskno, *elsub, *me, *color;
  int *newelsub;
  int nsub = expectedNumSubdomains;

  // Compute total profile size
  int totalSize = 0;
  for(i=0; i<decompNumSubdomains; ++i)
    totalSize += sizes[i];

  double percent = 0.95;
  int exact_num = etoe->numNonZeroP();
  int numele = etoe->csize();
  int numnod = ntoe->csize();
  elsub = new int [nsub+1];
  newelsub = new int[nsub+2];
  me = new int[exact_num+1];
  color = new int[numele];
  maskno = new int[numnod];

  //  initialize newelsub as a pointer into me 
  newelsub[0] = 0;

  int div = totalSize / nsub;
  int rem = totalSize % nsub;
  // fprintf(stderr,"Per Cluster %d Remainder %d\n",div,rem);

int div1 = exact_num / nsub;
int rem1 = exact_num % nsub;

  // Compute an a priori balanced distribution of element numbers
  for(i = 0; i < nsub; i++) {
    elsub[i] = (i < rem1) ? div1+1 : div1;
  }

  for(i = 0; i < nsub; i++) {
    newelsub[i+1] = newelsub[i] + elsub[i];
  }

  for(i = 0; i < nsub; i++) {
    elsub[i] = (i < rem) ? div+1 : div;
  }

  for(i = 0; i < numele; i++) color[i] = -1;

  for (i = 0; i < numnod; i++) maskno[i] = ntoe->num(i);

  ptr1 = 0;
  ptr2 = 0;
  int numsub = 0;
  ntot = 0;

  while(numsub < nsub && ntot <= percent*elsub[numsub]) {

    // Locate a starting (interface if possible) node
    // The part of code is in O(Numnod). It could be done in
    // O(I) (I = interface size) by handling properly the list of
    // interface nodes.

    minw = minw2 = 32000000;
    node = -1;
    for(i = 0; i < numnod; i++) {
      if(maskno[i] == 0) continue;
      if(maskno[i] == ntoe->num(i)) {
        if(maskno[i] < minw2) {
          minw2 = maskno[i];
          node2 = i;
        }
      } 
      else {
        if(maskno[i] < minw) {
          minw = maskno[i];
          node = i;
        }
      }
    }
 
    if(node < 0) node = node2;

    // Initialize the list of elements atached to the starting node
    ptr1 = ptr2;

    for(i = 0; i < ntoe->num(node) && ntot < percent*elsub[numsub]; i++) {
       // find an element attached to "node" 
       nelem = (*ntoe)[node][i];
       if(color[nelem] >= 0) continue;
       color[nelem] = numsub;
       me[ptr2] = nelem;
       ptr2++;
       ntot += sizes[nelem]; // MODIFICATION
       // reduce mask value for all nodes connected to this element 
       nomask(eton, maskno, nelem);
    }
 
    /* RECURSIVELY ADD TO LIST NEW ADJACENT ELEMENTS */

    while(ptr2 > ptr1 && ntot < percent*elsub[numsub]) {
      for(i = 0; i<etoe->num(me[ptr1]) && ntot<percent*elsub[numsub]; i++) {
        nelem = (*etoe)[me[ptr1]][i];
        if(color[nelem] >= 0) continue;
        color[nelem] = numsub;
        me[ptr2] = nelem;
        ntot += sizes[nelem]; // MODIFICATION
        ptr2++;
        nomask(eton, maskno, nelem);
      }
      ptr1++;
    }
 
    // MODIFICATION
    if(ntot >= percent*elsub[numsub] && ntot <= elsub[numsub]/percent) {
      numsub++;
      ntot = 0;
    } 
  }

  for(i = 0; i < numele; i++) { 
    if(etoe->num(i) && color[i] < 0) {
      me[ptr2++] = i; 
      color[i] = nsub;
    }
  }

  // delete [] maskno;
  // delete [] color;
  // delete [] elsub;

  Connectivity *cToS = new Connectivity(nsub, newelsub, me);

  return cToS;
}

bool 
ObjectOrdering::operator()(int n1, int n2)
{
  int nt1 = 0, nt2 = 0, t1, t2;
  do
    {
      while(nt1 < objToSub->num(n1) && subIsInClus[ (*objToSub)[n1][nt1] ] == false)
	nt1++;
      while(nt2 < objToSub->num(n2) && subIsInClus[ (*objToSub)[n2][nt2] ] == false)
	nt2++;
      t1 = ( nt1 < objToSub->num(n1) ) ? (*objToSub)[n1][nt1] : -1 ;
      t2 = ( nt2 < objToSub->num(n2) ) ? (*objToSub)[n2][nt2] : -1 ;
      if(t1 < t2)
	return(true);
      if(t2 < t1)
	return(false);
      nt1++;
      nt2++;
    } while (t1 >= 0);
  return( n1 < n2 );
}

struct DataEntry {
  TypeTag dataType;
  INT64BIT offset;
};

void Sower::write()
{
  for(int currentClusNum = 1 ; currentClusNum <= nCluster ; currentClusNum++)
    {
      // opening cluster file
      char filename[FILENAME_LENGTH];
      char clusterNumStr[10] = {'\0'};
      sprintf(clusterNumStr, "%u", currentClusNum);
      strcpy(filename, FILENAME_PREFIX);
      strcat(filename, clusterNumStr);
#ifdef SOWER_DEBUG
      cout << " ** Writing to file " << filename << endl;
#endif
      BinFileHandler file(filename, "wb");
      
      // counting number of data types in cluster
      int dataTypesNumber = 0;
      map<TypeTag,DataStruct*>::iterator it = entries.begin();
      for(; it!=entries.end(); it++) {
//cerr << "Passage" << (*it).second->getClusterToData(clusToSub)->num(currentClusNum-1) << endl;
//(*it).second->getClusterToData(clusToSub)->print();
	if((*it).second->getClusterToData(clusToSub)->num(currentClusNum-1) > 0) {
//cerr << "selection" << endl;
	  dataTypesNumber++;
        }
      }
      if(nSurfaces) dataTypesNumber++; //HB: make room for surfaces header/toc
#ifdef SOWER_DEBUG
      cout << " ** " << dataTypesNumber << " different data types will be written to this file." << endl;
#endif     
      
      // saving space for the table of content
      // Type of data ID / pointer to beginning of data / pointer to associated rangeset
      int tocSize = (sizeof(TypeTag)+sizeof(INT64BIT)+sizeof(INT64BIT)) * (dataTypesNumber+1);
      INT64BIT nextEntryOffset = file.tell() + tocSize;       /** offset in the data **/
      INT64BIT tocCurrentOffset = file.tell();  /** offset in the table of content **/
      INT64BIT curRangeSetLocation=0;
      // now writing each datatype
      it = entries.begin();
      for(; it!=entries.end(); it++)
	{
	  // oh no ! this datatype is not in this cluster
	  if((*it).second->getClusterToData(clusToSub)->num(currentClusNum-1) <= 0)
	    continue; 
	  // updating toc with entry for this datatype
	  file.seek(tocCurrentOffset);
             
          file.write(&(*it).first, 1);
	  file.write(&nextEntryOffset, 1);
	  tocCurrentOffset = file.tell();

	  // preparing to write all the data of this datatype
	  file.seek(nextEntryOffset);
	  (*it).second->write(currentClusNum, clusToSub, numSubdomains, file, curRangeSetLocation);
	  nextEntryOffset = file.tell();
	  
	  // adding location of offset of the Rangesets for this datatype to T O C
	  file.seek(tocCurrentOffset);
	  file.write(&curRangeSetLocation, 1);
	  tocCurrentOffset = file.tell();

	}

      //HB: create & write toc, range set & data for surfaces
      if(nSurfaces) { writeSurfaces(file, tocCurrentOffset, nextEntryOffset, curRangeSetLocation); } 
      
      // adding 0 entry to TOC to mark its end
      INT64BIT zero = 0;

      int zeroi = END_TYPE;
      file.seek(tocCurrentOffset);
      file.write(&zeroi,1);
      file.write(&zero,1);
      file.write(&zero,1);
      // cluster file closed by BinFileHandler destructor

    } // end cluster loop
}

size_t 
DataStruct::write(int clusNumber, Connectivity* clusToSub, int numSubdomains, 
                  BinFileHandler& file, INT64BIT& curRangeSetLocation)
{
  map<int, RangeSet*> rangeSet;  // where is each subdomain's data written ?

  // Connectivity * clusterToData = clusToSub->transcon(subToData);
  // Connectivity * clusterToData = getClusterToData(clusToSub); // PJSA
  if(!clusterToData) clusterToData = clusToSub->transcon(subToData);
//  cerr << "subToData = \n"; subToData->print();
//  cerr << "clusterToData = \n"; clusterToData->print();
  Connectivity * objToSub = subToData->reverse();
  bool * subIsInClus = new bool[numSubdomains];
  for(int i=0; i<numSubdomains; i++)
    subIsInClus[i] = false;
  for(int j = 0; j < clusToSub->num(clusNumber-1); j++) {
    subIsInClus[ ((*clusToSub)[clusNumber-1][j]) ] = true;
    rangeSet[ ((*clusToSub)[clusNumber-1][j]) ] = new RangeSet();
  }

  // sorting
  // the oject ids will be in an optimum order to be written
  objToSub->sortTargets();
  ObjectOrdering order(objToSub,subIsInClus);
  std::sort((*clusterToData)[clusNumber-1], (*clusterToData)[clusNumber-1] + clusterToData->num(clusNumber-1), order);
  // now let's write each object using its objectID
  std::list<int> currentSubs; // list of subdomains currently been written -- for use with rangeset
  int k;
  for(k = 0; k < clusterToData->num(clusNumber-1); k++) {
    int curObjID = (*clusterToData)[clusNumber-1][k];
    // range set
    std::list<int> nextCurrentSubs;
    int * listOfSubsObjIsIn = (*objToSub)[curObjID];
    int numOfSubsObjIsIn = objToSub->num(curObjID);
      
    for(int i = 0 ; i < numOfSubsObjIsIn; ++i) {
      if(subIsInClus[ listOfSubsObjIsIn[i] ]) { // sub is in cluster
        // std::list<int>::iterator it = find( currentSubs.begin(), currentSubs.end(), listOfSubsObjIsIn[i]); 
        std::list<int>::iterator it = std::find( currentSubs.begin(), currentSubs.end(), listOfSubsObjIsIn[i]); // PJSA: for sgi intel
        if(currentSubs.empty() || it == currentSubs.end()) { // this subdomain just appeared
          rangeSet[listOfSubsObjIsIn[i]]->start(file.tell(), curObjID);
          nextCurrentSubs.push_back(listOfSubsObjIsIn[i]); // this range might have to be ended at the next k
        }
        else { // this subdomain stays
          currentSubs.remove(listOfSubsObjIsIn[i]);
          nextCurrentSubs.push_back(listOfSubsObjIsIn[i]);
        }
      }
    }
    // all subdomains that are still in currentSubs have disapeared ! time to end them
    for(std::list<int>::iterator it = currentSubs.begin(); it != currentSubs.end(); ++it) {
      rangeSet[(*it)]->end((*clusterToData)[clusNumber-1][k-1]);
    } // can never be called for k = 0
    currentSubs.clear();
    currentSubs.merge(nextCurrentSubs);
    // end of range set
      
    // actual writing of data
    //file.write(&curObjID, 1);
    writeObjectData((*clusterToData)[clusNumber-1][k],file,curObjID);
  }
  // all subdomains that are still in currentSubs have disapeared ! time to end them
  for(std::list<int>::iterator it = currentSubs.begin(); it != currentSubs.end(); ++it)
    rangeSet[(*it)]->end((*clusterToData)[clusNumber-1][k-1]);
  curRangeSetLocation = file.tell();
  file.seek(file.tell() + sizeof(int));
  // writing range set for this datatype
  int realNumOfSubs = 0;
  for(map<int, RangeSet*>::iterator it = rangeSet.begin(); it != rangeSet.end(); ++it) {
    if(!(*it).second->empty()) {
      ++realNumOfSubs;
      int subNum = (*it).first;
      file.write(&subNum,1);
#ifdef SOWER_DEBUG
      (*it).second->print();
#endif
      (*it).second->write(file);
    }
  }
  INT64BIT tmp = file.tell();
  file.seek(curRangeSetLocation);
  file.write(&realNumOfSubs, 1);
  file.seek(tmp);
  return 0;  
}

#include <Driver.d/GeoSource.h>
size_t ElemsetIO::write(Elemset* eset, int index, BinFileHandler& file, int curObjID)
    {
      if ((*eset)[index])
	{//test JF for gap in elemset
	  int* nodes = (*eset)[index]->allNodes();
	  int numNodes = (*eset)[index]->numNodes();
	  int elType = (*eset)[index]->getElementType(); // PJSA
	  int globalID = ((*geoSource).localToGlobalElementsNumber)[curObjID];          // translation to global coords
#ifdef SOWER_DEBUG
	  std::cerr << "writing element index = " << index << " curObjID = " << curObjID << " globalID = " << globalID << std::endl;
#endif
	  file.write(&globalID, 1);
	  file.write(&elType, 1);
	  file.write(&numNodes, 1);
	  file.write(nodes, numNodes);
          double pressure = (*eset)[index]->getPressure();
          file.write(&pressure, 1);
          double preload = (*eset)[index]->getPreLoad();
          file.write(&preload, 1);
	}
      else 
	{
	  std::cerr << "Sower.h, void element in eset. index = " << index << std::endl;
	  int elType = -1;
	  file.write(&elType, 1);
	}
      return 0;
    }

//HB
void Sower::writeSurfaces(BinFileHandler& file, INT64BIT& tocCurrentOffset, INT64BIT& nextEntryOffset, INT64BIT& curRangeSetLocation)
{
#ifdef SOWER_DEBUG_SURFS
   cerr<<" ** Write surfaces"<<endl;
#endif
   // updating toc with entry for this datatype
   file.seek(tocCurrentOffset);
   int surfTagType = SURFACES;
   
   file.write(&surfTagType, 1); // i.e. DataType
   file.write(&nextEntryOffset, 1);
   tocCurrentOffset = file.tell();

   // preparing to write all the data of this datatype
   file.seek(nextEntryOffset);
   map<int, RangeSet*> surfRangeSet;  
   for(int isurf=0; isurf<nSurfaces; isurf++) { //assume the surface are "packed" in the Surfaces array
#ifdef SOWER_DEBUG_SURFS
   cerr<<" ** Write surface "<< isurf<<endl;
#endif
   surfRangeSet[isurf] = new RangeSet();
   surfRangeSet[isurf]->start(file.tell(), isurf);
   file.write(&isurf, 1);
   (*Surfaces)[isurf]->WriteSower(file);
   surfRangeSet[isurf]->end(isurf);
 }
  curRangeSetLocation = file.tell();
  file.seek(file.tell() + sizeof(int));
  // writing range set for this datatype
  int realNumOfSurfs = 0;
  for(map<int, RangeSet*>::iterator it = surfRangeSet.begin(); it != surfRangeSet.end(); ++it) {
    if(!(*it).second->empty()) {
      ++realNumOfSurfs;
      int surfNum = (*it).first;
//#ifdef SOWER_DEBUG_SURFS
//cerr << "Driver.d/Sower.C, surfNum = " << surfNum << endl;//-JFD
//#endif
      file.write(&surfNum,1);
#ifdef SOWER_DEBUG_SURFS
      (*it).second->print();
#endif
      (*it).second->write(file);
    }
  }
  INT64BIT tmp = file.tell();
  file.seek(curRangeSetLocation);
  file.write(&realNumOfSurfs, 1);
  file.seek(tmp);

  nextEntryOffset = file.tell();
  // adding location of offset of the Rangesets for this datatype to T O C
  file.seek(tocCurrentOffset);
  file.write(&curRangeSetLocation, 1);
  tocCurrentOffset = file.tell();
}
