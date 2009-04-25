#include <stdio.h>
#include <stdlib.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/Mpc.h>
#include <Utils.d/MFTT.h>
#include <Utils.d/BinFileHandler.h>
#include <Utils.d/Memory.h>
#include <Driver.d/SubDomainFactory.h>

#ifndef WINDOWS

#include <dlfcn.h>
#endif
#include <map>

#include <Driver.d/Access.h>

#ifdef SOWER_SURFS
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#endif

template<class Scalar>
void GeoSource::distributeBCs(GenSubDomain<Scalar> *&tmpSub, int *cl2LocNodeMap,
                                int *gl2ClusNodeMap)  
{
  MatrixTimers &mt = domain->getTimers();
  startTimerMemory(mt.distributeBCs,  mt.memoryDistBC);

  // get bc's for this subdomain
  BCond *subBC;
  BCond *subTextBC;

  int numLocDirichlet = getBC(dbc, numDirichlet, cl2LocNodeMap, subBC);

  // distribute dbc from textfile
  if (numTextDirichlet)  {
    int numLocTextDBC = getBC(textDBC, numTextDirichlet, cl2LocNodeMap,
                                subTextBC, gl2ClusNodeMap);
    if (numLocTextDBC)
      augmentBC(numLocTextDBC, subTextBC, subBC, numLocDirichlet);
  }

  if (numLocDirichlet > 0)
    tmpSub->setDirichlet(numLocDirichlet, subBC);

  int numLocNeuman = getBC(nbc, numNeuman, cl2LocNodeMap, subBC);

  if (numTextNeuman)  {
    int numLocTextNBC = getBC(textDBC, numTextNeuman, cl2LocNodeMap,
                                subTextBC, gl2ClusNodeMap);

    if (numLocTextNBC)
      augmentBC(numLocTextNBC, subTextBC, subBC, numLocNeuman);
  }

  if (numLocNeuman > 0)
    tmpSub->setNeuman(numLocNeuman, subBC);

  int numLocConvBC = getBC(cvbc, numConvBC, cl2LocNodeMap, subBC);
  if (numLocConvBC > 0)
    tmpSub->setConvBC(numLocConvBC, subBC);
 
  int numLocRadBC = getBC(rdbc, numRadBC, cl2LocNodeMap, subBC);
  if (numLocRadBC > 0)
    tmpSub->setRadBC(numLocRadBC, subBC);

  int numLocIDis = getBC(iDis, numIDis, cl2LocNodeMap, subBC);
  if (numLocIDis > 0)
    tmpSub->setIDis(numLocIDis, subBC);

  int numLocIDis6 = getBC(iDis6, numIDis6, cl2LocNodeMap, subBC);
  if (numLocIDis6 > 0)
    tmpSub->setIDis6(numLocIDis6, subBC);

  int numLocIVel = getBC(iVel, numIVel, cl2LocNodeMap, subBC);
  if (numLocIVel > 0)
    tmpSub->setIVel(numLocIVel, subBC);

  int numLocITemp = getBC(iTemp, numITemp, cl2LocNodeMap, subBC);
  if (numLocITemp > 0)
    tmpSub->setITemp(numLocITemp, subBC);

  //fprintf(stderr,"There are %d dbc's in sub %d \n", numLocDirichlet, tmpSub->subNum());

  stopTimerMemory(mt.distributeBCs,  mt.memoryDistBC);
}

#ifdef DISTRIBUTED
template<class Scalar>
void
GeoSource::distributeOutputNodes(GenSubDomain<Scalar> *&tmpSub,
                                 int *gl2ClNodeMap, int *cl2LocNodeMap)
{
  // count number of output nodes in this subdomain
  int iNode;
  int numLocOutNodes = 0;
  int *locOutputNodes = 0;
  int *locOutIndex = 0;

  for(iNode = 0; iNode < numNodalOutput; iNode++)  {
    if(outputNodes[iNode] <= maxGlobNode)  {
      int clusterNum = gl2ClNodeMap[outputNodes[iNode]];
      if(clusterNum < maxClusNode && clusterNum >= 0)
        if(cl2LocNodeMap[clusterNum] >= 0)
          numLocOutNodes++;
    }
  }

  if (numLocOutNodes)  {

    // allocate arrays to number of data in this subdomain
    locOutputNodes = new int[numLocOutNodes];
    locOutIndex = new int[numLocOutNodes];

    // populate arrays
    int count = 0;
    for (iNode = 0; iNode < numNodalOutput; iNode++)  {

      if (outputNodes[iNode] <= maxGlobNode)  {
        int clusterNum = gl2ClNodeMap[outputNodes[iNode]];
        if (clusterNum < maxClusNode && clusterNum >= 0)
          if (cl2LocNodeMap[clusterNum] >= 0)  {

            locOutputNodes[count] = outputNodes[iNode];
            locOutIndex[count] = outNodeIndex[iNode];

            // renumber to local node number
            locOutputNodes[count] = cl2LocNodeMap[clusterNum];

            // increment count
            count++;
          }
      }
    }
  }

  tmpSub->setOutputNodes(numLocOutNodes, locOutputNodes, locOutIndex);

}

//#define SOWER_DEBUG

template<class Scalar>
void
GeoSource::distributeOutputNodesX(GenSubDomain<Scalar> *tmpSub, Connectivity *nodeToSub)
{
  // count number of output nodes in this subdomain
  int iNode;
  int numLocOutNodes = 0;
  int *locOutputNodes = 0;
  int *locOutIndex = 0;

  for(iNode = 0; iNode < numNodalOutput; iNode++)  {
    if((*nodeToSub)[outputNodes[iNode]][0] == tmpSub->subNum())  {
      numLocOutNodes++;
    }
  }

  if (numLocOutNodes)  {
    // allocate arrays to number of data in this subdomain
    locOutputNodes = new int[numLocOutNodes];
    locOutIndex = new int[numLocOutNodes];
    // populate arrays
    int count = 0;
    for (iNode = 0; iNode < numNodalOutput; iNode++)  {
      if((*nodeToSub)[outputNodes[iNode]][0] == tmpSub->subNum())  {
        locOutputNodes[count] = tmpSub->globalToLocal(outputNodes[iNode]);
        // locOutputNodes[count] = outputNodes[iNode];
        locOutIndex[count] = outNodeIndex[iNode];
        count++;
      }
    }
  }

  tmpSub->setOutputNodes(numLocOutNodes, locOutputNodes, locOutIndex);
}
#endif

#ifndef SALINAS
#include <Driver.d/Sower.h>

template<class Scalar>
GenSubDomain<Scalar> *
GeoSource::readDistributedInputFiles(int localSubNum, int subNum)
{
#ifdef SOWER_DEBUG
  std::cerr<< "READING distributed input for subDomain " << subNum << std::endl;
#endif
  Sower sower;
 
  BinFileHandler *f = sower.openBinaryFile(subNum);

  // nodes
#ifdef SOWER_DEBUG
  std::cerr << std::endl << "--Reading Nodes" << std::endl;
#endif
  int *glNodeNums = 0;  // PJSA: this is the localToGlobal node numbering map, allocated in Sower::read
  CoordSet* cs = sower.template read<NodesIO>(*f, subNum, glNodeNums, true);
#ifdef SOWER_DEBUG
  if(cs) {
    for(int i = 0; i < cs->size(); i++)
      if((*cs)[i] != 0)
        std::cerr << "local " << i << ", global " << glNodeNums[i] << ", coords: " 
                  << (*cs)[i]->x << "," << (*cs)[i]->y << "," << (*cs)[i]->z << std::endl;
  }
#endif

  // elements
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Elements" << endl;
#endif
  int *glElemNums = 0; // PJSA: this is the localToGlobal element numbering map, allocated in Sower::read
  Elemset *ese = sower.template read<ElemsetIO>(*f, subNum, glElemNums);
#ifdef SOWER_DEBUG
  if(ese) ese->list();
#endif
  // pressure and preload
  prsflg = 1; prlflg = 1; 

  // materials
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Materials" << endl;
#endif
  int *glNums = 0;
  std::pair<int, map<int,StructProp>* > *rat0 = sower.template read<MatIO>(*f, subNum, glNums);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(rat0) {
    for(int l=0; l < rat0->first; ++l) {
      map<int,StructProp> *tm = rat0->second;
      cerr << l << ":" <<(*tm)[l].E << ","
           << (*tm)[l].A
           << "," << (*tm)[l].nu
           << "," << (*tm)[l].rho
           << "," << (*tm)[l].eh
           << "," << (*tm)[l].Ixx
           << "," << (*tm)[l].Iyy
           << "," << (*tm)[l].Izz
           << "," << (*tm)[l].c
           << "," << (*tm)[l].k
           << "," << (*tm)[l].Q
           << "," << (*tm)[l].W
           << "," << (*tm)[l].P
           << "," << (*tm)[l].Ta
           << "," << (*tm)[l].ymin
           << "," << (*tm)[l].ymax
           << "," << (*tm)[l].zmin
           << "," << (*tm)[l].zmax
           << (*tm)[l].kappaHelm << endl;
      }
  }
#endif

  // attributes
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Attributes" << endl;
#endif
  std::pair<int, map<int,Attrib>* > *rat = sower.template read<AttribIO>(*f, subNum, glNums);
  // composites
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Composites" << endl;
#endif
  /*std::pair<int,map<int,CoefData*>* > *rcd =*/ sower.template read<CompositeCIO>(*f, subNum, glNums);
  /*std::pair<int,map<int,LayInfo*>* > *rli =*/ sower.template read<CompositeLIO>(*f, subNum, glNums);
  /*std::pair<int,map<int,double*>* > *rcf =*/ sower.template read<CFramesIO>(*f, subNum, glNums);

  if(rat) {
    for(int iAttr=0; iAttr<rat->first; ++iAttr) {
      int locElemNum = (*rat->second)[iAttr].nele;
      int att = (*rat->second)[iAttr].attr;
      StructProp *prop = &(*rat0->second)[att];
      (*ese)[locElemNum]->setProp(prop);
      int cmp_att = (*rat->second)[iAttr].cmp_attr;
      if(cmp_att > -1) { }
      // PJSA: warning need special treatment for composite and phantom elements here  
    }
  }
  if(glNums) { delete [] glNums; glNums = 0; } // not used

  GenSubDomainFactory<Scalar>* pFactory;
  pFactory = GenSubDomainFactory< Scalar >::getFactory();
  GenSubDomain<Scalar> *subd = pFactory->
    createSubDomain(*domain, localSubNum, cs, ese, glNodeNums, glElemNums, subNum);
  // don't delete glNodeNums and glElemNums, they now belong to the SubDomain object

  // neuman bounday conditions
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Neuman Boundary Conditions (forces)" << std::endl;
#endif
  list<BCond *> *fo = sower.template read<BCDataIO<FORCES_TYPE> >(*f, subNum, glNums);
  if(fo) subd->setNeumanBC(fo);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(fo) {
    for(list<BCond *>::iterator it = fo->begin(); it != fo->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  //boffset
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Beam offsets" << endl;
#endif
  vector<OffsetData> *offsets2 = sower.template read<BoffsetIO>(*f, subNum, glNums);
  if(offsets2) {
    for(vector<OffsetData>::iterator it = (*offsets2).begin() ; it!= (*offsets2).end() ;++it)
      {
        offsets.push_back(*it);
      }
  }

  //eframes
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Eframes" << endl;
#endif
  vector<EFrameData> *efd2 = sower.template read<EFrameIO>(*f,subNum,glNums);
  if(efd2) {
    for(vector<EFrameData>::iterator it = (*efd2).begin();it!= (*efd2).end();++it)
      {
        double fr[9];
        fr[0] = (*it).frame[0][0];
        fr[1] = (*it).frame[0][1];
        fr[2] = (*it).frame[0][2];
        fr[3] = (*it).frame[1][0];
        fr[4] = (*it).frame[1][1];
        fr[5] = (*it).frame[1][2];
        fr[6] = (*it).frame[2][0];
        fr[7] = (*it).frame[2][1];
        fr[8] = (*it).frame[2][2];
        geoSource->setFrame((*it).elnum, fr);
      }
  }

  // dirichlet boundary conditions
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Dirichlet Boundary=disp" << endl;
#endif
  list<BCond *> *fo2 = sower.template read<BCDataIO<DISPLACEMENTS_TYPE> >(*f, subNum, glNums);
  if(fo2) subd->setDirichletBC(fo2);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(fo2) {
    for(list<BCond *>::iterator it = fo2->begin(); it != fo2->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // convection
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Convection" << endl;
#endif
  list<BCond *> *clist = sower.template read<BCDataIO<CONV_TYPE> >(*f, subNum, glNums);
  if(clist) subd->setConvection(clist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(clist) {
    for(list<BCond *>::iterator it = clist->begin(); it != clist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // radiation
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Radiation" << endl;
#endif
  list<BCond *> *rlist = sower.template read<BCDataIO<RAD_TYPE> >(*f, subNum, glNums);
  if(rlist) subd->setRadiation(rlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(clist) {
    for(list<BCond *>::iterator it = rlist->begin(); it != rlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // initial displacements
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Initial Displacement" << endl;
#endif
  list<BCond *> *idlist = sower.template read<BCDataIO<IDISP_TYPE> >(*f, subNum, glNums);
  if(idlist) subd->setInitialDisplacement(idlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(idlist) {
    for(list<BCond *>::iterator it = idlist->begin(); it != idlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // initial displacements (6 column)
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Initial Displacement 6" << endl;
#endif
  list<BCond *> *id6list = sower.template read<BCDataIO<IDISP6_TYPE> >(*f, subNum, glNums);
  if(id6list) subd->setInitialDisplacement6(id6list);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(id6list) {
    for(list<BCond *>::iterator it = id6list->begin(); it != id6list->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // initial velocities
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Initial Velocity" << endl;
#endif
  list<BCond *> *ivlist = sower.template read<BCDataIO<IVEL_TYPE> >(*f, subNum, glNums);
  if(ivlist) subd->setInitialVelocity(ivlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(ivlist) {
    for(list<BCond *>::iterator it = ivlist->begin(); it != ivlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // initial temperatures
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Initial Temperature" << endl;
#endif
  list<BCond *> *itlist = sower.template read<BCDataIO<ITEMP_TYPE> >(*f, subNum, glNums);
  if(itlist) subd->setInitialTemperature(itlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(itlist) {
    for(list<BCond *>::iterator it = itlist->begin(); it != itlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  if(claw) subd->setClaw(claw->fileName, claw->routineName);

  // claw->sensors
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Sensors (claw)" << endl;
#endif
  list<BCond *> *sensorlist = sower.template read<BCDataIO<SENSOR_TYPE> >(*f, subNum, glNums);
  if(sensorlist) subd->setSensor(sensorlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(sensorlist) {
    for(list<BCond *>::iterator it = sensorlist->begin(); it != sensorlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif

  // claw->actuator
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Actuator (claw)" << endl;
#endif
  list<BCond *> *actuatorlist = sower.template read<BCDataIO<ACTUATOR_TYPE> >(*f, subNum, glNums);
  if(actuatorlist) subd->setActuator(actuatorlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(actuatorlist) {
    for(list<BCond *>::iterator it = actuatorlist->begin(); it != actuatorlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif
                                                                                                 
  // claw->userDisp
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Usdd (claw)" << endl;
#endif
  list<BCond *> *usddlist = sower.template read<BCDataIO<USDD_TYPE> >(*f, subNum, glNums);
  if(usddlist) subd->setUsdd(usddlist); 
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(usddlist) {
    for(list<BCond *>::iterator it = usddlist->begin(); it != usddlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif
  
  // claw->userForce
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Usdf (claw)" << endl;
#endif
  list<BCond *> *usdflist = sower.template read<BCDataIO<USDF_TYPE> >(*f, subNum, glNums);
  if(usdflist) subd->setUsdf(usdflist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(usdflist) {
    for(list<BCond *>::iterator it = usdflist->begin(); it != usdflist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->val << endl;
    }
  }
#endif
  if (claw) claw->makeGlobalClaw(subd->getClaw());//Warning, claw for Geosource is "empty"

  // complex dirichlet
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Complex Dirichlet" << endl;
#endif
  list<ComplexBCond *> *cdlist = sower.template read<ComplexBCDataIO<HDIR_TYPE> >(*f, subNum, glNums);
  if(cdlist) subd->setComplexDirichletBC(cdlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(cdlist) {
    for(list<ComplexBCond *>::iterator it = cdlist->begin(); it != cdlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->reval 
                << "," << (*it)->imval << endl;
    }
  }
#endif

  // complex neuman
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Complex Neuman" << endl;
#endif
  list<ComplexBCond *> *cnlist = sower.template read<ComplexBCDataIO<HNEU_TYPE> >(*f, subNum, glNums);
  if(cnlist) subd->setComplexNeumanBC(cnlist);
  if(glNums) { delete [] glNums; glNums = 0; } // not used
#ifdef SOWER_DEBUG
  if(cnlist) {
    for(list<ComplexBCond *>::iterator it = cnlist->begin(); it != cnlist->end(); ++it) {
      std::cerr << (*it)->nnum << " : " << (*it)->dofnum << "," << (*it)->reval 
                << "," << (*it)->imval << endl;
    }
  }
#endif

  //ATDDNB & HDNB
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading DNB (ATDDNB or HDNB)" << endl;
#endif
  list<SommerElement *> *dnblist = sower.template read<SommerDataIO<DNB_TYPE> >(*f, subNum, glNums);
  if(dnblist) subd->setDnb(dnblist);
  if(glNums) { delete[] glNums; glNums = 0; } //not used
#ifdef SOWER_DEBUG
  if(dnblist) {
    int i, numN = 0;
    for(list<SommerElement *>::iterator it = dnblist->begin(); it != dnblist->end(); ++it) {
      //get type and list of nodes
      numN = (*it)->numNodes();
      std::cerr << " node: ";
      for (i = 0 ; i < numN-1 ; i++)
        std::cerr << (*it)->getNode(i) << ", ";
      std::cerr << (*it)->getNode(numN-1) << endl;
    }
  }
#endif

  //LMPC
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading LMPC" << endl;
#endif
  std::pair<int,ResizeArray<LMPCons*>* > *mpc = sower.template read<LMPCIO >(*f, subNum, glNums);
  if(mpc) {
    for(int i=0; i<mpc->first; ++i) domain->addLMPC((*mpc->second)[i]);
  }
  if(glNums) { delete[] glNums; glNums = 0; }
#ifdef SOWER_DEBUG
  if(mpc) {
    cerr << "numLMPC = " << mpc->first << endl;
    for(int i=0; i<mpc->first; ++i) {
      (*mpc->second)[i]->print();
    }
  }
#endif
  
  //TETT
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading TETT" << endl;
#endif
  /*std::pair<int, ResizeArray<MFTTData*>* >* tett =*/ sower.template read<MFTTDataIO<TETT_TYPE> >(*f, subNum, glNums);

  //YMTT
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading YMTT" << endl;
#endif
  /*std::pair<int, ResizeArray<MFTTData*>* >* ymtt =*/ sower.template read<MFTTDataIO<YMTT_TYPE> >(*f, subNum, glNums);

  //ATDROB & HSCB
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading Scatter (ATDROB & HSCB)" << endl;
#endif
  list<SommerElement *> *scatlist = sower.template read<SommerDataIO<SCAT_TYPE> >(*f, subNum, glNums);
  if(scatlist) subd->setScat(scatlist);
  if(glNums) {delete[] glNums; glNums = 0; } //not used
#ifdef SOWER_DEBUG
  if(scatlist) {
    int i, numN = 0;
    for(list<SommerElement *>::iterator it = scatlist->begin(); it != scatlist->end(); ++it) {
      //get type and list of nodes
      numN = (*it)->numNodes();
      std::cerr << " node: ";
      for (i = 0 ; i < numN-1 ; i++)
        std::cerr << (*it)->getNode(i) << ", ";
      std::cerr << (*it)->getNode(numN-1) << endl;
    }
  }
#endif

  //ATDARB & HARB
#ifdef SOWER_DEBUG
  cerr << endl << "--Reading ARB (ATDARB & HARB)" << endl;
#endif
  list<SommerElement *> *arblist = sower.template read<SommerDataIO<ARB_TYPE> >(*f, subNum, glNums);
  if(arblist)  subd->setArb(arblist);
  if(glNums) {delete[] glNums; glNums = 0; } //not used
#ifdef SOWER_DEBUG
  if(arblist) {
    int i, numN = 0;
    for(list<SommerElement *>::iterator it = arblist->begin(); it != arblist->end(); ++it) {
      //get type and list of nodes
      numN = (*it)->numNodes();
      std::cerr << " node: ";
      for (i = 0 ; i < numN-1 ; i++)
        std::cerr << (*it)->getNode(i) << ", ";
      std::cerr << (*it)->getNode(numN-1) << endl;
    }
  }
#endif

  cleanUp();
  return subd;
}
#else
template<class Scalar>
GenSubDomain<Scalar> *
GeoSource::readDistributedInputFiles(int localSubNum, int subNum)
{
    return NULL;    //CRW
}
#endif

