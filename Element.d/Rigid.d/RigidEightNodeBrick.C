#include <Element.d/Rigid.d/RigidEightNodeBrick.h>
#include <Element.d/Joint.d/ConstantDistanceConstraint.h>

RigidEightNodeBrick::RigidEightNodeBrick(int *_nn)
{
  nSubElems = 18;
  subElems = new Element * [nSubElems];
  int indices[18][2] = {{0,1},
                        {1,2},
                        {2,3},
                        {4,5},
                        {5,6},
                        {6,7},
                        {0,2},
                        {4,6},
                        {0,3},
                        {4,7},
                        {0,4},
                        {1,5},
                        {2,6},
                        {3,7},
                        {0,5},
                        {1,6},
                        {2,7},
                        {0,7}};
  for(int i = 0; i < nSubElems; ++i) {
    subElems[i] = new ConstantDistanceConstraint(indices[i]);
  }
  initialize(8, _nn);
}

