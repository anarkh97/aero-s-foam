#include <Driver.d/GeoSource.h>

#include <Utils.d/BinaryOutputFile.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <cstddef>
#include <stdexcept>
#include <cassert>

// Private functions

void GeoSource::getOutputFileName(char *result, int fileId, int clusterId, int iter) {
  const char *suffix = computeClusterSuffix(clusterId + 1, clusToSub->csize());
  if (outLimit > -1) {
    sprintf(result, "%s%d_%s", oinfo[fileId].filename, iter/outLimit, suffix);
  } else {
    sprintf(result, "%s%s", oinfo[fileId].filename, suffix);
  }
  delete[] suffix;
}

inline
int
GeoSource::getHeaderNameBytes(int fileId) const {
  return headLen[fileId] * sizeof(char);
}

BinFileHandler *
GeoSource::openBinaryOutputFile(int fileId, int clusterId, int iter, const char *flag)
{
  char outfileName[64];
  getOutputFileName(outfileName, fileId, clusterId, iter);
  return new BinFileHandler(outfileName, flag);
}

void GeoSource::outputHeader(BinFileHandler &file, int dim, int fileId)
{
  char headDescrip[200];
  const int headerLength = getHeaderDescriptionAndLength(headDescrip, fileId);
  const std::string headerDesc(headDescrip, headerLength);

  const int dataType = oinfo[fileId].dataType;
  if (dataType != 1 && dataType != 2) {
    throw std::logic_error("Bad Data Type");
  }

  const int numClusItems = (dataType == 1) ? numClusNodes : numClusElems;
  
  writeHeaderToBinaryOutputFile(file, dataType, headerDesc, numClusItems, dim);
}

void 
GeoSource::writeArrayToBinFile(const double *data, int dataSize, int subId, int inDataOffset, int fileId, 
                               int iterRank, int resultRank, double timeStamp, int inStateDataCount, int clusterItemCount) {
  const int clusterId = (*subToClus)[subId][0];
  const char *appendFlag = "ws+";
  std::auto_ptr<BinFileHandler> binFile(openBinaryOutputFile(fileId, clusterId, iterRank, appendFlag));
  
  const int firstSubInCluster = (*clusToSub)[clusterId][0];
  const bool doSerialPart = (firstSubInCluster == subId);

  writePartialResultToBinaryOutputFile(*binFile, data, dataSize, inDataOffset, doSerialPart, 
                                       resultRank, timeStamp, inStateDataCount, clusterItemCount,
                                       getHeaderNameBytes(fileId));
}

// Public functions

void GeoSource::createBinaryOutputFile(int fileId, int glSub, int iter)
{
   if(!oinfo[fileId].PodRomfile) {
  // Open file for first time and write header
  if(oinfo[fileId].interval != 0) {
    if (binaryOutput) {
      // Determine which cluster subdomain is in
      const int clusterId = (*subToClus)[glSub][0];
      const char *truncateFlag = "w";
      std::auto_ptr<BinFileHandler> binFile(openBinaryOutputFile(fileId, clusterId, iter, truncateFlag));

      // Write header information
      outputHeader(*binFile, oinfo[fileId].dim, fileId);
    } else { // ASCII output
      std::ofstream outfile;
      outfile.open(oinfo[fileId].filename, std::ofstream::out | std::ofstream::trunc);
      char headDescrip[200];
      const int headerLength = getHeaderDescriptionAndLength(headDescrip, fileId);
      headDescrip[headerLength] = '\0'; // Ensure we have a valid C-style string
      outfile << headDescrip;
      outfile.close();
    }
  }
 }
}

void GeoSource::outputRange(int fileId, int *globalIndex, int nData, int glSub, int offset, int iter)
{
  if(!oinfo[fileId].PodRomfile) {
   if (binaryOutput) {
    const int clusterId = (*subToClus)[glSub][0];
    const char *appendFlag = "ws+";
    std::auto_ptr<BinFileHandler> file(openBinaryOutputFile(fileId, clusterId, iter, appendFlag));

    const int dataType = oinfo[fileId].dataType;
    if (dataType != 1 && dataType != 2) {
      throw std::logic_error("Bad Data Type");
    }

    if (dataType == 1) {
      // Nodal data
      const bool compressedNumbering = (domain->outFlag == 1);
      
      std::vector<int> compressedIndex;
      for (int iNode = 0; iNode < nData; ++iNode) {
        const int nodeIndex = globalIndex[iNode];
        const bool isMasterNode = (nodeIndex >= 0);

        if (isMasterNode && nodeIndex < numNodes /*nodes.size()*/) { // note: nodes.size() is not available when using "binaryinput on"
          compressedIndex.push_back(compressedNumbering ? (domain->nodeTable[nodeIndex] - 1) : nodeIndex);
        }
      }

      const int* compressedIndexArray = compressedIndex.empty() ? NULL : &compressedIndex[0];
      writePartialRangeToBinaryOutputFile(*file, compressedIndexArray, compressedIndex.size(), offset, getHeaderNameBytes(fileId)); 
    } else {
      // Elemental data
      writePartialRangeToBinaryOutputFile(*file, globalIndex, nData, offset, getHeaderNameBytes(fileId)); 
    }
  }
 }
}


void 
GeoSource::writeNodeScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber,
                                 int iter, int numRes, double time, int numComponents, int *glNodeNums)
{
  const int group = oinfo[fileNumber].groupNumber;
  if (binaryOutput && group == -1) { // Group output is always ASCII
    writeArrayToBinFile(data, numData, glSub, offset, fileNumber, iter, numRes, time, numComponents * numClusNodes, numClusNodes);
  } else {
    std::ofstream outfile;
    outfile.open(oinfo[fileNumber].filename, std::ios_base::in | std::ios_base::out);
    outfile.precision(oinfo[fileNumber].precision);

    int numComponentsPlus = (group == -1) ? numComponents : numComponents + 4; // group output: allow for NODENUMBER, X0, Y0, Z0
    int numNodesPlus = (group == -1) ? (domain->outFlag ? domain->exactNumNodes : numNodes) : nodeGroup[group].size();

    long timeOffset = long(headLen[fileNumber])  // header including std::endl
                      + long(numRes-1)*(3 + oinfo[fileNumber].width + 1)  // 3 spaces + time(s) + std::endl for all previous timesteps
                      + long(numRes-1)*long(numNodesPlus)*(numComponentsPlus*(2+oinfo[fileNumber].width)+1);
                                                                     // 2 spaces + first_component + ... + 2 spaces + last_component + std::endl
                                                                     // for each node for all previous timesteps

    outfile.setf(std::ios_base::showpoint | std::ios_base::right | std::ios_base::scientific | std::ios_base::uppercase);
    int glNode_prev = -1;
    if(glSub == 0) { // if first subdomain in cluster, write time
      outfile.seekp(timeOffset);
      outfile.width(3+oinfo[fileNumber].width);
      outfile << time << std::endl;
      // fix for gaps in node numbering (note: this isn't required for group output or when domain->outFlag != 0)
      // the first subdomain writes zeros for all unasigned nodes
      // TODO for "binaryinput on", Domain->nodeToElem is NULL this fix doesn't work: we need to precompute a list
      // of unassigned nodes, somehow, without using nodeToElem.
      if(group == -1 && domain->outFlag == 0 && domain->getNodeToElem() != NULL) {
        int counter = 0;
        for(int i=0; i<numNodes /*nodes.size()*/; ++i) { // note: nodes.size() is not available when using "binaryinput on"
          if(domain->getNodeToElem()->num(i) == 0) {
            int glNode = i;
            if(glNode-glNode_prev != 1) { // need to seek in file for correct position to write next node
              long relativeOffset = long(glNode-glNode_prev-1)*(numComponents*(2+oinfo[fileNumber].width) + 1);
              outfile.seekp(relativeOffset, std::ios_base::cur);
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
    timeOffset += (3+oinfo[fileNumber].width+1); // 3 spaces + width + 1 for std::endl
    if(glSub != 0) outfile.seekp(timeOffset);
  
    int k = 0;
    for(int i = 0; i < numData/numComponents; ++i) {
      while(true) { if(glNodeNums[k] == -1) k++; else break; }
      if(glNodeNums[k] >= numNodes /*nodes.size()*/) { k++; continue; } // don't print "internal" nodes eg for rigid beams
                                                               // note: nodes.size() is not available when using "binaryinput on"
      int glNode = (domain->outFlag) ? domain->nodeTable[glNodeNums[k]]-1 : glNodeNums[k]; k++;
      if(group != -1) {
        std::list<int>::iterator it = nodeGroup[group].begin();
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
        long relativeOffset = long(glNode-glNode_prev-1)*(numComponentsPlus*(2+oinfo[fileNumber].width) + 1);
        outfile.seekp(relativeOffset, std::ios_base::cur);
      }
      if(group != -1) { // print NODENUMBER, X0, Y0, Z0
        outfile.width(2+oinfo[fileNumber].width);
        outfile << glNodeNums[k-1]+1;
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

  delete[] data;
}

void 
GeoSource::writeElemScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber,
                                 int iter, int numRes, double time, int totData, int *glElemNums)  
{
  if (binaryOutput) {
    writeArrayToBinFile(data, numData, glSub, offset, fileNumber, iter, numRes, time, totData, numClusElems);
  } else {
    std::ofstream outfile;
    outfile.open(oinfo[fileNumber].filename, std::ios_base::in | std::ios_base::out);
    outfile.precision(oinfo[fileNumber].precision);

    long timeOffset = headLen[fileNumber]  // header including std::endl
                      + long(numRes-1)*(3 + oinfo[fileNumber].width + 1) 
                      + long(numRes-1)*(long(totData)*(2+oinfo[fileNumber].width)+nElem);

    outfile.setf(std::ios_base::showpoint | std::ios_base::right | std::ios_base::scientific | std::ios_base::uppercase);
    if(glSub == 0) { // if first subdomain in cluster, write time
      outfile.seekp(timeOffset);
      outfile.width(3+oinfo[fileNumber].width);
      outfile << time << std::endl;
    }
    timeOffset += (3+oinfo[fileNumber].width+1); // 3 spaces + width + 1 for std::endl

    int numNodesPerElem = totData/nElem;
    if(glSub==0 && time == 0) std::cerr << " *** WARNING: text file output for element stresses in distributed mode is only correctly\n"
                                        << " ***          implemented for models in which all elements have the same number of nodes \n";
     
    int k = 0;
    for(int i=0; i<numData/numNodesPerElem; ++i) {
      int glElem = glElemNums[i]; 
      long totalOffset = timeOffset + long(glElem)*((2+oinfo[fileNumber].width)*numNodesPerElem + 1);
      outfile.seekp(totalOffset);
      for(int j=0; j<numNodesPerElem; ++j) { outfile.width(2+oinfo[fileNumber].width); outfile << data[k++]; }
      outfile << std::endl;
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

  delete[] data;
}

