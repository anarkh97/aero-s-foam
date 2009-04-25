#ifndef _NLHEXAHEDRAL_H_
#define _NLHEXAHEDRAL_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/3DShapeFunction.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/StrainDispEvaluator.h>



class NLHexahedral : public GaussIntgElement {

    int n[8];

    Material *material;
    
    bool linearKinematics;


protected:
   int getNumGaussPoints();
   void getGaussPointAndWeight(int i, double *point, double &weight);
   ShapeFunction *getShapeFunction();
   StrainEvaluator *getStrainEvaluator();
   Material *getMaterial();

public:
   NLHexahedral(int *nd, bool isLinKin);
   int numNodes() { return 8; }
   int numDofs() { return 24; }
   PrioInfo examine(int sub, MultiFront *); // dec
   void renum(int *);
   void   markDofs(DofSetArray &);
   int*   dofs(DofSetArray &, int *p=0);
   int*   nodes(int * = 0);
   void updateStates(Node *nodes, double *states, double *un, double *unp) {}
   void setProp(StructProp *);
   void setMaterial(Material *);
   int getTopNumber();
};

#endif
