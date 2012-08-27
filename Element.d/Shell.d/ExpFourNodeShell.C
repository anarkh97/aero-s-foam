#include <Element.d/Beam.d/EulerBeam.h>
#include <Element.d/Shell.d/ExpFourNodeShell.h>
#include <Corotational.d/utilities.h>

ExpFourNodeShell::ExpFourNodeShell(int *nodenums)
{
  int i,j,k;
  nn = new int[4];
  for(i=0; i<4; ++i) nn[i] = nodenums[i];
 
  nSubElems = 6;
  subElems = new Element * [6];
  subElemNodes = new int * [6];
  subElemDofs = new int * [6];

  subElemNodes[0] = new int[2];
  subElemNodes[0][0] = 0; subElemNodes[0][1] = 1;

  subElemNodes[1] = new int[2];
  subElemNodes[1][0] = 1; subElemNodes[1][1] = 2;

  subElemNodes[2] = new int[2];
  subElemNodes[2][0] = 2; subElemNodes[2][1] = 3; 

  subElemNodes[3] = new int[2];
  subElemNodes[3][0] = 3; subElemNodes[3][1] = 0; 

  subElemNodes[4] = new int[2];
  subElemNodes[4][0] = 0; subElemNodes[4][1] = 2; 

  subElemNodes[5] = new int[2];
  subElemNodes[5][0] = 1; subElemNodes[5][1] = 3; 


  for(i=0; i<6; ++i) {
    int tmp[2];
    subElemDofs[i] = new int[12];
    for(j=0; j<2; ++j) {
      int nij = subElemNodes[i][j];
      tmp[j] = nodenums[nij]; // global node numbers
      for(k=0;k<6;++k) {
        subElemDofs[i][6*j+k] = 6*nij+k;   
      }
    }
    subElems[i] = new EulerBeam(tmp);
    subElems[i]->setGlNum(-1);
  }
  nnodes = 4;
  ndofs = 24;
}

void 
ExpFourNodeShell::buildFrame(CoordSet& cs) 
{
  double a[3], b[3], normal[3];
  a[0] = cs.getNode(nn[1]).x - cs.getNode(nn[0]).x;
  a[1] = cs.getNode(nn[1]).y - cs.getNode(nn[0]).y;
  a[2] = cs.getNode(nn[1]).z - cs.getNode(nn[0]).z;
  b[0] = cs.getNode(nn[3]).x - cs.getNode(nn[0]).x;
  b[1] = cs.getNode(nn[3]).y - cs.getNode(nn[0]).y;
  b[2] = cs.getNode(nn[3]).z - cs.getNode(nn[0]).z;
  crossprod(a, b, normal);
  normalize(normal);

  EFrame* elemframe = new EFrame[6];
  for(int i=0; i<6; ++i) {
    EFrame& theFrame = elemframe[i];
    theFrame[0][0] = cs.getNode(subElemNodes[i][1]).x - cs.getNode(subElemNodes[i][0]).x;
    theFrame[0][1] = cs.getNode(subElemNodes[i][1]).y - cs.getNode(subElemNodes[i][0]).y;
    theFrame[0][2] = cs.getNode(subElemNodes[i][1]).z - cs.getNode(subElemNodes[i][0]).z;
    theFrame[2][0] = normal[0];
    theFrame[2][1] = normal[1];
    theFrame[2][2] = normal[2];
    normalize(theFrame[0]);
    crossprod(theFrame[0],theFrame[2],theFrame[1]);
    normalize(theFrame[1]);
    subElems[i]->setFrame(&elemframe[i]);
  }
}

bool
ExpFourNodeShell::isShell()
{
  return true;
}

int
ExpFourNodeShell::getTopNumber()
{
  return 188;
}

