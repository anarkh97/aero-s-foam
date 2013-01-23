#include <Element.d/Element.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Element.d/Sommerfeld.d/Line2SommerBC.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/Triangle6SommerBC.h>
#include <HelmAxi.d/LineAxiSommer.h>
#include <HelmAxi.d/Line2AxiSommer.h>
#include <HelmAxi.d/ScatterData.h>

ScatterData *globalScatter = 0;

ScatterData::ScatterData() : subNodes(0), scatEle(0) {

  numNode = 0;
  numElem = 0;
  numFFP  = 0;

}

void ScatterData::addScatter(SommerElement *ele) {

 scatEle[numElem++] = ele;

}

void
ScatterData::addScatterElem(int num, int etype, double sommerConst, int nnodes, 
                            int *n) {

   SommerElement *ele;

   switch(etype) {
     case 1:
       ele = new LineSommerBC(n[0], n[1]);
       addScatter(ele);
       break;
     case 3:
       ele = new TriangleSommerBC(n[0], n[1], n[2]);
       addScatter(ele);
       break;
     case 11:
       ele = new LineAxiSommer(n[0], n[1]);
       addScatter(ele);
       break;
     case 12:
       ele = new Line2AxiSommer(n[0], n[1], n[2]);
       addScatter(ele);
       break;
/*
     case 2:
       ele = new Line2SommerBC(n[0], n[1], n[2], sommerConst);
       addScatter(ele);
       break;
     case 4:
       ele = new QuadSommerBC(n[0], n[1], n[2], n[3], sommerConst);
       addScatter(ele);
       break;
     case 6:
       ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5], 
                                   sommerConst);
       addScatter(ele);
       break;
*/
     default:
       return;
   }

}

void ScatterData::subtractNode(int n) {

 subNodes[numNode++] = n;

}

void ScatterData::setFFP(int n) {

 numFFP = n;

}


