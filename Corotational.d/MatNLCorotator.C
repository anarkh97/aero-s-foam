#include <Corotational.d/MatNLCorotator.h>
#include <Element.d/NLElement.h>
#include <Corotational.d/GeomState.h>

MatNLCorotator::MatNLCorotator(MatNLElement *_ele, bool _own)
  : ele(_ele), own(_own)
{

}

MatNLCorotator::~MatNLCorotator()
{
  if(own) delete ele;
}

void
MatNLCorotator::getStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &C0,
                                 FullSquareMatrix &elk, double *f, double dt, double t)
{
  int *nn = new int[ele->numNodes()];
  ele->nodes(nn);
  Node *nodes = new Node[ele->numNodes()];
  for(int i = 0; i < ele->numNodes(); ++i) nodes[i] = *(C0[nn[i]]);

  double *dispn = new double[ele->numDofs()];
  double *staten;
  if(refState) {
    for(int i = 0; i < ele->numNodes(); ++i) {
      dispn[3*i+0] = (*refState)[nn[i]].x-nodes[i].x;
      dispn[3*i+1] = (*refState)[nn[i]].y-nodes[i].y;
      dispn[3*i+2] = (*refState)[nn[i]].z-nodes[i].z;
    }
    if(refState == &curState) {
      double *staten_ = refState->getElemState(ele->getGlNum());
      staten = new double[ele->numStates()];
      for(int i = 0; i < ele->numStates(); ++i) staten[i] = staten_[i];
    }
    else
      staten = refState->getElemState(ele->getGlNum());
  }
  else {
    for(int i = 0; i < ele->numDofs(); ++i) dispn[i] = 0.0;
    staten = new double[ele->numStates()];
    for(int i = 0; i < ele->numStates(); ++i) staten[i] = 0.0;
  }

  double *dispnp = new double[ele->numDofs()];
  for(int i = 0; i < ele->numNodes(); ++i) {
    dispnp[3*i+0] = curState[nn[i]].x-nodes[i].x;
    dispnp[3*i+1] = curState[nn[i]].y-nodes[i].y;
    dispnp[3*i+2] = curState[nn[i]].z-nodes[i].z;
  }
  double *statenp = curState.getElemState(ele->getGlNum());

  ele->integrate(nodes, dispn, staten, dispnp, statenp, elk, f, dt);
  for(int i = 0; i < ele->numDofs(); ++i) f[i] = -f[i];

  delete [] nn;
  delete [] nodes;
  delete [] dispn;
  if(!refState || refState == &curState) delete [] staten;
  delete [] dispnp;
}

void
MatNLCorotator::getInternalForce(GeomState *refState, GeomState &curState, CoordSet &C0,
                                 FullSquareMatrix &, double *f, double dt, double t)
{
  int *nn = new int[ele->numNodes()];
  ele->nodes(nn);
  Node *nodes = new Node[ele->numNodes()];
  for(int i = 0; i < ele->numNodes(); ++i) nodes[i] = *(C0[nn[i]]);

  double *dispn = new double[ele->numDofs()];
  double *staten;
  if(refState) {
    for(int i = 0; i < ele->numNodes(); ++i) {
      dispn[3*i+0] = (*refState)[nn[i]].x-nodes[i].x;
      dispn[3*i+1] = (*refState)[nn[i]].y-nodes[i].y;
      dispn[3*i+2] = (*refState)[nn[i]].z-nodes[i].z;
    }
    if(refState == &curState) {
      double *staten_ = refState->getElemState(ele->getGlNum());
      staten = new double[ele->numStates()];
      for(int i = 0; i < ele->numStates(); ++i) staten[i] = staten_[i];
    }
    else 
      staten = refState->getElemState(ele->getGlNum());
  }
  else {
    for(int i = 0; i < ele->numDofs(); ++i) dispn[i] = 0.0;
    staten = new double[ele->numStates()];
    for(int i = 0; i < ele->numStates(); ++i) staten[i] = 0.0;
  }

  double *dispnp = new double[ele->numDofs()];
  for(int i = 0; i < ele->numNodes(); ++i) {
    dispnp[3*i+0] = curState[nn[i]].x-nodes[i].x;
    dispnp[3*i+1] = curState[nn[i]].y-nodes[i].y;
    dispnp[3*i+2] = curState[nn[i]].z-nodes[i].z;
  }
  double *statenp = curState.getElemState(ele->getGlNum());

  ele->integrate(nodes, dispn, staten, dispnp, statenp, f, dt);
  for(int i = 0; i < ele->numDofs(); ++i) f[i] = -f[i];

  delete [] nn;
  delete [] nodes;
  delete [] dispn;
  if(!refState || refState == &curState) delete [] staten;
  delete [] dispnp;
}

void
MatNLCorotator::extractDeformations(GeomState &curState, CoordSet &C0, double *u, int &nlflag)
{
  int *nn = new int[ele->numNodes()];
  ele->nodes(nn);
  for(int i = 0; i < ele->numNodes(); i++)  {
    u[3*i+0] = curState[nn[i]].x - C0[nn[i]]->x;
    u[3*i+1] = curState[nn[i]].y - C0[nn[i]]->y;
    u[3*i+2] = curState[nn[i]].z - C0[nn[i]]->z;
  }
  delete [] nn;
  nlflag = 2;
}

void
MatNLCorotator::getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                              GeomState *refState, CoordSet &C0, int strIndex, 
                              int surface, double *ndTemps, double ylayer, 
                              double zlayer, int avgnum, int measure)
{
  int *nn = new int[ele->numNodes()];
  ele->nodes(nn);
  Node *nodes = new Node[ele->numNodes()];
  for(int i = 0; i < ele->numNodes(); ++i) nodes[i] = *(C0[nn[i]]);

  double *dispn = new double[ele->numDofs()];
  double *staten;
  if(refState) {
    for(int i = 0; i < ele->numNodes(); ++i) {
      dispn[3*i+0] = (*refState)[nn[i]].x-nodes[i].x;
      dispn[3*i+1] = (*refState)[nn[i]].y-nodes[i].y;
      dispn[3*i+2] = (*refState)[nn[i]].z-nodes[i].z;
    }
    staten = refState->getElemState(ele->getGlNum());
  }
  else {
    for(int i = 0; i < ele->numDofs(); ++i) dispn[i] = 0.0;
    staten = new double[ele->numStates()];
    for(int i = 0; i < ele->numStates(); ++i) staten[i] = 0.0;
  }

  double *dispnp = new double[ele->numDofs()];
  for(int i = 0; i < ele->numNodes(); ++i) {
    dispnp[3*i+0] = curState[nn[i]].x-nodes[i].x;
    dispnp[3*i+1] = curState[nn[i]].y-nodes[i].y;
    dispnp[3*i+2] = curState[nn[i]].z-nodes[i].z;
  }
  double *statenp = curState.getElemState(ele->getGlNum());
  
  int indexMap[6] = { 0, 4, 8, 1, 5, 2 };
  int numPoints = (avgnum == -1) ? ele->getNumGaussPoints() : ele->numNodes();
  switch(strIndex) {
    case 0 : case 1 : case 2 : case 3 : case 4 : case 5 : { // SXX=0,SYY=1,SZZ=2,SXY=3,SYZ=4,SXZ=5
      double (*result)[9] = new double[numPoints][9];
      double *statenp_tmp = new double[ele->numStates()];
      ele->getStressTens(nodes, dispn, staten, dispnp, statenp_tmp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i][indexMap[strIndex]];
      }
      delete [] result;
      delete [] statenp_tmp;
    } break;
    case 6 : { // VON
      double *result = new double[numPoints];
      double *statenp_tmp = new double[ele->numStates()];
      ele->getVonMisesStress(nodes, dispn, staten, dispnp, statenp_tmp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i];
      }
      delete [] result;
      delete [] statenp_tmp;
    } break;
    case 7 : case 8 : case 9 : case 10 : case 11 : case 12 : { // EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12
      double (*result)[9] = new double[numPoints][9];
      ele->getStrainTens(nodes, dispnp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i][indexMap[strIndex-7]];
      }
      if(strIndex == 10 || strIndex == 11 || strIndex == 12) {
        for(int i = 0; i < numPoints; ++i) {
          stress[i] *= 2; // convert to "engineering strain"
        }
      }
      delete [] result;
    } break;
    case 13 : { // STRAINVON
      double *result = new double[numPoints];
      ele->getVonMisesStrain(nodes, dispnp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i];
      }
      delete [] result;
    } break;
    case 18 : { // EQPLSTRN
      double *result = new double[numPoints];
      ele->getEquivPlasticStrain(statenp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) stress[i] = result[i];
      delete [] result;
    } break;
    case 19 : case 20 : case 21 : case 22 : case 23 : case 24 : { // BACKSXX=19,BACKSYY=20,BACKSZZ=21,
                                                                  // BACKSXY=22,BACKSYZ=23,BACKSXZ=24
      double (*result)[9] = new double[numPoints][9];
      ele->getBackStressTens(statenp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i][indexMap[strIndex-19]];
      }
      delete [] result;
    } break;
    case 25 : case 26 : case 27 : case 28 : case 29 : case 30 : { // PLASTICEXX=25,PLASTICEYY=26,PLASTICEZZ=27,
                                                                  // PLASTICEXY=28,PLASTICEYZ=29,PLASTICEXZ=30
      double (*result)[9] = new double[numPoints][9];
      ele->getPlasticStrainTens(statenp, result, avgnum);
      for(int i = 0; i < numPoints; ++i) {
        stress[i] = result[i][indexMap[strIndex-25]];
      }
      delete [] result;
    } break;
  }

  weight = 1;

  delete [] nn;
  delete [] nodes;
  delete [] dispn;
  if(!refState) delete [] staten;
  delete [] dispnp;
}

void
MatNLCorotator::updateStates(GeomState *refState, GeomState &curState, CoordSet &C0)
{
  int *nn = new int[ele->numNodes()];
  ele->nodes(nn);
  Node *nodes = new Node[ele->numNodes()];
  for(int i = 0; i < ele->numNodes(); ++i) nodes[i] = *(C0[nn[i]]);

  double *dispn = new double[ele->numDofs()];
  if(refState) {
    for(int i = 0; i < ele->numNodes(); ++i) {
      dispn[3*i+0] = (*refState)[nn[i]].x-nodes[i].x;
      dispn[3*i+1] = (*refState)[nn[i]].y-nodes[i].y;
      dispn[3*i+2] = (*refState)[nn[i]].z-nodes[i].z;
    }
  }
  else {
    for(int i = 0; i < ele->numDofs(); ++i) dispn[i] = 0.0;
  }

  double *dispnp = new double[ele->numDofs()];
  for(int i = 0; i < ele->numNodes(); ++i) {
    dispnp[3*i+0] = curState[nn[i]].x-nodes[i].x;
    dispnp[3*i+1] = curState[nn[i]].y-nodes[i].y;
    dispnp[3*i+2] = curState[nn[i]].z-nodes[i].z;
  }
  double *state = curState.getElemState(ele->getGlNum());
  //for(int i=0; i<103; ++i) std::cerr << state[i] << " "; std::cerr << std::endl;
  
  ele->updateStates(nodes, state, dispn, dispnp);

  delete [] nn;
  delete [] nodes;
  delete [] dispn;
  delete [] dispnp;
}

int
MatNLCorotator::getNumGaussPoints()
{
  return ele->getNumGaussPoints();
}
