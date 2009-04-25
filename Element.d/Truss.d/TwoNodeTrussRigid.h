#ifndef _TWONODETRUSSRIGID_H_
#define _TWONODETRUSSRIGID_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>

class TwoNodeTrussRigid : public Element , public Corotator{

        int nn[3];
public:

	TwoNodeTrussRigid(int*);

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&,double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);


        void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        int              numInternalNodes() { return 1; }
        void             setInternalNodes(int *in) { nn[2] = in[0]; }

        Corotator *getCorotator(CoordSet &, double *, int, int)
           { return this; }

    void getStiffAndForce(GeomState &, CoordSet &,
                                  FullSquareMatrix &, double *);

    // ONLY FOR BUCKLING (EIGEN)
    // ONLY NONLINEAR TERM OF K
     void formGeometricStiffness(GeomState &, CoordSet &,
                                 FullSquareMatrix &, double *) { };

     double* getOriginalStiffness() { return 0; };

     void extractDeformations(GeomState &geomState, CoordSet &cs,
                              double *vld, int &nlflag) { };

     void getNLVonMises(Vector&, Vector& weight,
                        GeomState &, CoordSet &, int) { };

     void getNLAllStress(FullM&, Vector&,
                         GeomState &, CoordSet &, int) { };

     double getElementEnergy(GeomState &, CoordSet &) { return 0.0; };

     void extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                        double *vlr) { };

     int  getTopNumber() { return 101; }
     bool isRigidElement() { return true; }
};
#endif
