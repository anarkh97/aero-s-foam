#ifndef _MPCELEMENT_H_
#define _MPCELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
class LMPCons;

class MpcElement : public Element , public Corotator
{
    int nnodes;  // number of nodes
    int *nn;  // node numbers
    int ndofs; // number of dofs
    LMPCons *mpc;
    int *renumTable;

  public:
    MpcElement(LMPCons *mpc);
    virtual ~MpcElement()
      { if(nn) delete [] nn;  if(renumTable) delete [] renumTable; };

    void setProp(StructProp *p) { };
    bool isMpcElement() { return true; }

    void renum(int *table);

    FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix imagStiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int  numDofs() { return ndofs; }
    int  numNodes() { return nnodes; }
    int* nodes(int * = 0);
    int  numInternalNodes() { return 1; }
    void setInternalNodes(int *in) { nn[nnodes] = in[0]; nnodes++; }

    Corotator *getCorotator(CoordSet &, double *, int, int)
       { return this; }

    void getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &, double *);

     void formGeometricStiffness(GeomState &, CoordSet &, FullSquareMatrix &, double *)
       { cerr << "formGeometricStiffness() not implemented for MpcElement \n"; };

     double* getOriginalStiffness() { return 0; };

     void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag)
       { /*cerr << "extractDeformations() not implemented for MpcElement \n";*/ };

     void getNLVonMises(Vector&, Vector& weight, GeomState &, CoordSet &, int)
       { /*cerr << "getNLVonMises() not implemented for MpcElement \n";*/ };

     void getNLAllStress(FullM&, Vector&, GeomState &, CoordSet &, int)
       { /*cerr << "getNLAllStress() not implemented for MpcElement \n";*/ };

     double getElementEnergy(GeomState &, CoordSet &) { return 0.0; };

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs, double *vlr)
       { /*cerr << "extractRigidBodyMotion() not implemented for MpcElement \n";*/ };

     int  getTopNumber() { return 501; }
     int  numTopNodes() { return nnodes-1; }

     PrioInfo examine(int sub, MultiFront *mf);

     bool isSafe() { return false; }
};
#endif

