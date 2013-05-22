#ifndef _FSIELEMENT_H_
#define _FSIELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
class LMPCons;

class FsiElement : public Element
{
    int nnodes;  // number of nodes
    int *nn;  // node numbers
    int ndofs; // number of dofs
    int *renumTable;
    LMPCons *fsi;

  public:
    FsiElement(LMPCons *mpc);
    virtual ~FsiElement()
      { if(nn) delete [] nn;  if(renumTable) delete [] renumTable; };

    void setProp(StructProp *p) { };
    bool isFsiElement() { return true; }
// JLchange: only valid for the case of one-one connection fsi 
    int fsiFluidNode()  { return nn[nnodes-1]; }
    int fsiStrutNode()  { return nn[0]; } 

    LMPCons* cons() { return fsi; }

    void renum(int *table);
        void renum(EleRenumMap&);

    FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix imagStiffness(CoordSet&, double *kel, int flg=1);
    FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);

    void markDofs(DofSetArray &);
    int* dofs(DofSetArray &, int *p=0);
    int  numDofs() { return ndofs; }
    int  numNodes() { return nnodes; }
    int* nodes(int * = 0);
    int  numInternalNodes() { return 0; }

    int  getTopNumber() { return 502; }
    int  numTopNodes() { return nnodes; }

    PrioInfo examine(int sub, MultiFront *mf);
// JLchange    bool isSafe() { return false; }
};
#endif

