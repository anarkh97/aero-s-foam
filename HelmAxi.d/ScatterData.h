#ifndef _SCATTERDATA_H_
#define _SCATTERDATA_H_

#include <Utils.d/resize_array.h>

class CoordSet;
class SommerElement;

class ScatterData {

public : 

  int numNode;
  int numElem;

  int numFFP;

  ResizeArray<int> subNodes;
  ResizeArray<SommerElement *> scatEle;   

  ScatterData();
  void addScatter(SommerElement *ele);
  void addScatterElem(int, int, double, int, int *);
  void subtractNode(int n);
  void setFFP(int);

};

extern ScatterData *globalScatter;

#endif
