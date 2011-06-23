#include <Utils.d/BinFileHandler.h>
#include <Driver.d/GeoSource.h>
#include <fstream>

#define TEXT_OUTPUT

void GeoSource::createBinaryOutputFile(int iInfo, int glSub, int iter)
{
  // open file for 1st time and write header
  if(oinfo[iInfo].interval != 0) {
    if(binaryOutput) {
      // determine which cluster subdomain is in
      int clusNum = (*subToClus)[glSub][0];

      // determine output file name suffixes
      char outfileName[32];
      char *suffix = computeClusterSuffix(clusNum+1, clusToSub->csize());
      if(outLimit > -1) sprintf(outfileName, "%s%d_%s", oinfo[iInfo].filename, iter/outLimit, suffix); 
      else sprintf(outfileName, "%s%s", oinfo[iInfo].filename, suffix);

      delete [] suffix;

      // open binary output file for the first time
      BinFileHandler binFile(outfileName, "ws");

      // write header information
      outputHeader(binFile, oinfo[iInfo].dim, iInfo);
    } 
    else {
      ofstream outfile;
      outfile.open(oinfo[iInfo].filename, ofstream::out | ofstream::trunc);
      if(headLen == 0) headLen = new int[numOutInfo];
      char headDescrip[200];   // description of data
      getHeaderDescription(headDescrip, iInfo);
      headLen[iInfo] = strlen(headDescrip);
      for(int i=0; i< headLen[iInfo]; ++i) outfile << headDescrip[i];
      outfile.close();
    }
  }
}

void GeoSource::outputHeader(BinFileHandler &file, int dim, int fileNum)
{
  // write data description
  if(headLen == 0) headLen = new int[numOutInfo];
  char headDescrip[200];   // description of data
  getHeaderDescription(headDescrip, fileNum);
  headLen[fileNum] = strlen(headDescrip);

  int dataType = oinfo[fileNum].dataType;
  if(dataType != 1 && dataType != 2)
    fprintf(stderr," *** ERROR: Bad Data Type: %d\n", dataType);

  file.write(&dataType, 1);
  file.write(headLen+fileNum, 1);
  file.write(headDescrip, headLen[fileNum]);

  // write number of cluster data
  if(dataType == 1)   // nodal data
    file.write(&numClusNodes, 1);
  else  // elemental data
    file.write(&numClusElems, 1);

  // write dimension of data
  file.write(&dim, 1);

  // save spot for number of results
  file.write(&dim, 1);
}

void GeoSource::setHeaderLen(int fileNum)
{
  // write data description
  if(headLen == 0) headLen = new int[numOutInfo];
  char headDescrip[200];   // description of data
  getHeaderDescription(headDescrip, fileNum);
  headLen[fileNum] = strlen(headDescrip);
}

void GeoSource::outputRange(int fileNum, int *flag, int nData, int glSub, int offset, int iter)
{
  if(binaryOutput) {
    int clusNum = (*subToClus)[glSub][0];
    BinFileHandler *file = openBinaryOutputFile(clusNum, fileNum, iter);

    // skip according to offset
    file->seek( (offset+5) * sizeof(int) + headLen[fileNum]*sizeof(char) );

    int dataType = oinfo[fileNum].dataType;
    if(dataType != 1 && dataType != 2)
      fprintf(stderr," *** ERROR: Bad Data Type: %d\n", dataType);

    if(dataType == 1)  {
      int *flag2 = (domain->outFlag) ? new int[nData] : flag;
      int range = 0;
      int rStart = 0;
      for(int iNode = 0; iNode < nData; iNode++)  {
        if(domain->outFlag) flag2[iNode] = (flag[iNode] >=0) ? domain->nodeTable[flag[iNode]]-1 : flag[iNode];
        if(flag[iNode] >= 0 && flag[iNode] < nodes.size()) range++; // now skipping internal nodes
        else {
          if(range > 0) {
            file->write(flag2+rStart, range);
            range = 0;
          }
          rStart = iNode+1;
        }
      }
      if(range > 0) file->write(flag2+rStart, range);
      if(domain->outFlag) delete [] flag2;
    }
    else  {
      file->write(flag, nData);
    }
    delete file;
  }
}

BinFileHandler *
GeoSource::openBinaryOutputFile(int clusNum, int fileNumber, int iter)
{
  // determine output file name suffixes
  char outfileName[32];
  char *suffix = computeClusterSuffix(clusNum+1, clusToSub->csize());
  if(outLimit > -1) sprintf(outfileName, "%s%d_%s", oinfo[fileNumber].filename, iter/outLimit, suffix); // PJSA 3-31-06
  else sprintf(outfileName, "%s%s", oinfo[fileNumber].filename, suffix);

  delete [] suffix;

  // open binary output file
  const char *flag = "ws+";
  return new BinFileHandler(outfileName, flag);
}

void 
GeoSource::writeNodeScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber, 
                                 int iter, int numRes, double time, int numComponents, int *glNodeNums)
{
  int group = oinfo[fileNumber].groupNumber;
  if(binaryOutput && group == -1) { // group output is always ascii
    // get number of cluster files to write to
    int clusNum = (*subToClus)[glSub][0];

    // open output file
    BinFileHandler *binFile = openBinaryOutputFile(clusNum, fileNumber, iter);

    BinFileHandler::OffType timeOffset = headLen[fileNumber]*sizeof(char) 
                                         + (5+numClusNodes)*sizeof(int)
                                         + (numRes-1)*(numComponents*numClusNodes+1)*sizeof(double);

    if((*clusToSub)[clusNum][0] == glSub) {  // if first subdomain in cluster 
      // update number of results
      binFile->seek(4*sizeof(int)+headLen[fileNumber]*sizeof(char));
      binFile->write(&numRes, 1);
      // write time tag
      binFile->seek(timeOffset);
      binFile->write(&time, 1);
    }

    // seek to the correct position in file
    BinFileHandler::OffType totalOffset = timeOffset + 1*sizeof(double) + offset*sizeof(double);
    binFile->seek(totalOffset);

    binFile->write(data, numData);

    delete binFile;
  }
  else {
    ofstream outfile;
    outfile.open(oinfo[fileNumber].filename,  ios_base::in | ios_base::out);
    outfile.precision(oinfo[fileNumber].precision);

    int numComponentsPlus = (group == -1) ? numComponents : numComponents + 4; // group output: allow for NODENUMBER, X0, Y0, Z0
    int numNodesPlus = (group == -1) ? (domain->outFlag ? domain->exactNumNodes : numNodes) : nodeGroup[group].size();

    long timeOffset = headLen[fileNumber]  // header including endl
                      + (numRes-1)*(3 + oinfo[fileNumber].width + 1)  // 3 spaces + time(s) + endl for all previous timesteps
                      + (numRes-1)*numNodesPlus*(numComponentsPlus*(2+oinfo[fileNumber].width)+1);
                                                                     // 2 spaces + first_component + ... + 2 spaces + last_component + endl
                                                                     // for each node for all previous timesteps

    outfile.setf(ios_base::showpoint | ios_base::right | ios_base::scientific | ios_base::uppercase);
    int glNode_prev = -1;
    if(glSub == 0) { // if first subdomain in cluster, write time
      outfile.seekp(timeOffset);
      outfile.width(3+oinfo[fileNumber].width);
      outfile << time << endl;
      // fix for gaps in node numbering (note: this isn't required for group output or when domain->outFlag != 0)
      // the first subdomain writes zeros for all unasigned nodes
      if(group == -1) {
        int counter = 0;
        for(int i=0; i<nodes.size(); ++i) {
          if(domain->getNodeToElem()->num(i) == 0 && domain->outFlag == 0) {
            int glNode = i;
            if(glNode-glNode_prev != 1) { // need to seek in file for correct position to write next node
              long relativeOffset = (glNode-glNode_prev-1)*(numComponents*(2+oinfo[fileNumber].width) + 1);
              outfile.seekp(relativeOffset, ios_base::cur);
            }
            for(int j=0; j<numComponents; ++j) { 
              outfile.width(2+oinfo[fileNumber].width);
              outfile << double(0);
            }
            outfile << "\n";
            glNode_prev = glNode;
            counter++;
          }
        }
      }
    }
    timeOffset += (3+oinfo[fileNumber].width+1); // 3 spaces + width + 1 for endl
    if(glSub != 0) outfile.seekp(timeOffset);
  
    int k = 0;
    for(int i = 0; i < numData/numComponents; ++i) {
      while(true) { if(glNodeNums[k] == -1) k++; else break; }
      if(glNodeNums[k] >= nodes.size()) continue; // don't print "internal" nodes eg for rigid beams
      int glNode = (domain->outFlag) ? domain->nodeTable[glNodeNums[k]]-1 : glNodeNums[k]; k++;
      if(group != -1) {
        list<int>::iterator it = nodeGroup[group].begin();
        int grNode = 0;
        while(it != nodeGroup[group].end()) {
          int inode = *it;
          if(inode == glNode) break;
          it++; grNode++;
        }
        if(it == nodeGroup[group].end()) continue;
        else glNode = grNode;
      }
      if(glNode-glNode_prev != 1) { // need to seek in file for correct position to write next node
        //long totalOffset = timeOffset + glNode*(numComponents*(2+oinfo[fileNumber].width) + 1);
        //outfile.seekp(totalOffset);
        long relativeOffset = (glNode-glNode_prev-1)*(numComponentsPlus*(2+oinfo[fileNumber].width) + 1);
        outfile.seekp(relativeOffset, ios_base::cur);
      }
      if(group != -1) { // print NODENUMBER, X0, Y0, Z0
        outfile.width(2+oinfo[fileNumber].width);
        outfile << glNodeNums[k-1];
        outfile.width(2+oinfo[fileNumber].width);
        outfile << nodes[glNodeNums[k-1]]->x;
        outfile.width(2+oinfo[fileNumber].width);
        outfile << nodes[glNodeNums[k-1]]->y;
        outfile.width(2+oinfo[fileNumber].width);
        outfile << nodes[glNodeNums[k-1]]->z;
      }
      for(int j = 0; j < numComponents; ++j) { 
        outfile.width(2+oinfo[fileNumber].width); 
        outfile << data[i*numComponents+j];
      }
      outfile << "\n";
      glNode_prev = glNode;
    }
    outfile.close(); 
  }
}

void
GeoSource::writeNodeScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, 
                                 int iter, int numRes, double time, int numComponents, int *glNodeNums)
{
  double *data = new double[numData];
  switch(oinfo[fileNumber].complexouttype) {
    default:
    case OutputInfo::realimag :
      // extract real part of data and write to file
      for(int i=0; i<numData; ++i) data[i] = complexData[i].real();
      writeNodeScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes-1, time, numComponents, glNodeNums);
      // extract imaginary part of data and write to file
      for(int i=0; i<numData; ++i) data[i] = complexData[i].imag();
      writeNodeScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes, time, numComponents, glNodeNums);
      break;
    case OutputInfo::modulusphase :
      // compute moduli and write to file
      for(int i=0; i<numData; ++i) data[i] = abs(complexData[i]);
      writeNodeScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes-1, time, numComponents, glNodeNums);
      // compute phases and write to file
      for(int i=0; i<numData; ++i) data[i] = arg(complexData[i]);
      writeNodeScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes, time, numComponents, glNodeNums);
      break;
    case OutputInfo::animate :
      int nco = oinfo[fileNumber].ncomplexout;
      double phi = 0;
      double incr = 2.0*PI/double(nco);
      for(int j=0; j<nco; ++j) {
        for(int i=0; i<numData; ++i) data[i] = abs(complexData[i])*cos(arg(complexData[i])-phi);
        writeNodeScalarToFile(data, numData, glSub, offset, fileNumber, iter, nco*(numRes-1)+j+1, phi, numComponents, glNodeNums);
        phi += incr;
      }
      break;
  }

  delete [] data;
}

void 
GeoSource::writeElemScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber,
                                 int iter, int numRes, double time, int totData, int *glElemNums)  
{
  if(binaryOutput) {
    // get number of cluster files to write to
    int clusNum = (*subToClus)[glSub][0];

    // open output file
    BinFileHandler *binFile = openBinaryOutputFile(clusNum, fileNumber, iter);

    BinFileHandler::OffType timeOffset = headLen[fileNumber]*sizeof(char)
                                         + (5+numClusElems)*sizeof(int)
                                         + (numRes-1)*(totData+1)*sizeof(double);

    if((*clusToSub)[clusNum][0] == glSub) { // if first subdomain in cluster
      // update number of results
      binFile->seek(4*sizeof(int)+headLen[fileNumber]*sizeof(char));
      binFile->write(&numRes, 1);
      // write time tag
      binFile->seek(timeOffset);
      binFile->write(&time, 1);
    }

    // seek to data offset postion to write data
    BinFileHandler::OffType totalOffset = timeOffset + 1*sizeof(double) + offset*sizeof(double);
    binFile->seek(totalOffset);

    binFile->write(data, numData);
    delete binFile;
  }
  else {
    ofstream outfile;
    outfile.open(oinfo[fileNumber].filename,  ios_base::in | ios_base::out);
    outfile.precision(oinfo[fileNumber].precision);

    long timeOffset = headLen[fileNumber]  // header including endl
                      + (numRes-1)*(3 + oinfo[fileNumber].width + 1) 
                      + (numRes-1)*(totData*(2+oinfo[fileNumber].width)+nElem);

    outfile.setf(ios_base::showpoint | ios_base::right | ios_base::scientific | ios_base::uppercase);
    if(glSub == 0) { // if first subdomain in cluster, write time
      outfile.seekp(timeOffset);
      outfile.width(3+oinfo[fileNumber].width);
      outfile << time << endl;
    }
    timeOffset += (3+oinfo[fileNumber].width+1); // 3 spaces + width + 1 for endl

    int numNodesPerElem = totData/nElem;
    if(glSub==0) cerr << " *** WARNING: text file output for element stresses in distributed mode is only correctly\n"
                      << " ***          implemented for models in which all elements have the same number of nodes \n";
     
    int k = 0;
    for(int i=0; i<numData/numNodesPerElem; ++i) {
      int glElem = glElemNums[i]; 
      long totalOffset = timeOffset + glElem*((2+oinfo[fileNumber].width)*numNodesPerElem + 1);
      outfile.seekp(totalOffset);
      for(int j=0; j<numNodesPerElem; ++j) { outfile.width(2+oinfo[fileNumber].width); outfile << data[k++]; }
      outfile << endl;
    }
    outfile.close();
  }
}

void
GeoSource::writeElemScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, 
                                 int iter, int numRes, double time, int totData, int *glElemNums)
{
  double *data = new double[numData];

  switch(oinfo[fileNumber].complexouttype) {
    default:
    case OutputInfo::realimag :
      // extract real part of data and write to file
      for(int i=0; i<numData; ++i) data[i] = complexData[i].real();
      writeElemScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes-1, time, totData, glElemNums);
      // extract imaginary part of data and write to file
      for(int i=0; i<numData; ++i) data[i] = complexData[i].imag();
      writeElemScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes, time, totData, glElemNums);
      break;
    case OutputInfo::modulusphase :
      // compute moduli and write to file
      for(int i=0; i<numData; ++i) data[i] = abs(complexData[i]);
      writeElemScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes-1, time, totData, glElemNums);
      // compute phases and write to file
      for(int i=0; i<numData; ++i) data[i] = arg(complexData[i]);
      writeElemScalarToFile(data, numData, glSub, offset, fileNumber, iter, 2*numRes, time, totData, glElemNums);
      break;
    case OutputInfo::animate :
      int nco = oinfo[fileNumber].ncomplexout;
      double phi = 0;
      double incr = 2.0*PI/double(nco);
      for(int j=0; j<nco; ++j) {
        for(int i=0; i<numData; ++i) data[i] = abs(complexData[i])*cos(arg(complexData[i])-phi);
        writeElemScalarToFile(data, numData, glSub, offset, fileNumber, iter, nco*(numRes-1)+j+1, phi, totData, glElemNums);
        phi += incr;
      }
      break;
  }

  delete [] data;
}

