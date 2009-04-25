#ifndef _SOMMERELEMENT_H_
#define _SOMMERELEMENT_H_

#include <Element.d/Element.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/ComplexD.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/FullRectMatrix.h>
#include <Driver.d/Mpc.h>

class Domain;

class SommerElement : public Element {

private:
        SommerElement(const SommerElement &);

protected:
       static bool first; //HB for debugging coupled pb

public:
	Element *el;//Adjacent 3D element
	Element *el2;
        Domain *dom;
        complex<double> soundSpeed; // for wet scattering rhs
        SommerElement(Element* _el=0, Domain* _dom=0) { el = _el; dom = _dom; 
                                                        el2 = 0; } //HB
        virtual int numNodes()=0;
	virtual int getNode(int);
	virtual int* getNodes();
        virtual int* nodes(int * = 0);
        virtual int numDofs()=0;
        virtual int numWetDofs();
        virtual int numSolidDofs();
        virtual int dim();
        virtual int* dofs(DofSetArray &, int *p=0)=0;
        virtual int* wetDofs(DofSetArray &, int *p=0);
        virtual int* solidDofs(DofSetArray &, int *p=0);
        virtual void markDofs(DofSetArray &);
	virtual void renum(int *);
        virtual SommerElement* clone();

        virtual int nFaceCorners() { return 0; }
        virtual int* faceCorners() { return 0; }

        virtual void flipNormal();
        virtual void checkAndSetNormal(CoordSet& cs);
        virtual void center(CoordSet& cs, double* c);
        virtual void findBothEle (Connectivity *nodeToElem, int *eleTouch,
                             int *eleCount, int myNum, Elemset* eset, int *ie);
        virtual int findAndSetBothEle (CoordSet& cs, Elemset &eset,
             Connectivity *nodeToElem, int *eleTouch, int *eleCount, int myNum);

	virtual int findEle(Connectivity *nodeToEle, int *eleTouch, 
                            int *eleCount, int myNum, Elemset* eset=0,
                            int it = 0);
        virtual int findAndSetEle(CoordSet& cs,Elemset &eset,
             Connectivity *nodeToEle, int *eleTouch, int *eleCount, int myNum,
             int it = 0);

        virtual FullSquareMatrix sommerMatrix(CoordSet&);
        virtual FullSquareMatrix turkelMatrix(CoordSet&);

        virtual FullSquareMatrix refinedSommerMatrix(CoordSet&);
        //virtual FullSquareMatrix surfStiffMatrix(CoordSet&);
        virtual FullSquareMatrix HSommerMatrix(CoordSet&);
        //virtual FullSquareMatrix HKSommerMatrix(CoordSet&);
        virtual FullSquareMatrix interfMatrixConsistent(CoordSet&);
        virtual FullSquareMatrix interfMatrixLumped(CoordSet&);
        virtual FullSquareMatrix sommerMatrix(CoordSet&, double *)=0;
        virtual GenStackFSFullMatrix<double> wetInterfaceMatrix(CoordSet&,
                                                                double*);
        virtual void wetInterfaceLMPC(CoordSet &cs, LMPCons *lmpc, int nd);
        virtual void ffp(CoordSet &cs, int numFFP, double *dirFFP,
                         complex<double> *sol, complex<double> *ffpv);

        virtual FullRectMatrix transferMatrix(CoordSet&, double *);
        virtual FullSquareMatrix turkelMatrix(CoordSet&, double *);
        virtual FullSquareMatrix refinedSommerMatrix(CoordSet&, double *);
        //virtual FullSquareMatrix surfStiffMatrix(CoordSet&, double *);
        virtual FullSquareMatrix HSommerMatrix(CoordSet&, double *);
        //virtual FullSquareMatrix HKSommerMatrix(CoordSet&, double *);
        virtual FullSquareMatrix interfMatrixConsistent(CoordSet&, double *);
        virtual FullSquareMatrix interfMatrixLumped(CoordSet&, double *);
	virtual void BT2(CoordSet& cs, double *e, double *f, double *g,
                     double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d);
        virtual void BT2n(CoordSet& cs, double *e, double *f, double *g,
                     double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n);
        virtual void sphereBT2(CoordSet& cs, double r, double k, ComplexD *d);
        virtual void ellipsoidBT2(CoordSet& cs, double a, double b, double k, ComplexD *d);
        virtual void sommerMatrixEllipsoid(CoordSet &cs, double kappa, double H[3] , double K[3], ComplexD *d);

        virtual void neumVector(CoordSet&,Vector&, int pflag=0);
        virtual void neumVector(CoordSet&,ComplexVector&,
                                double,double,double,double, int pflag=0);
// RT: obsolete? and I need it for something else
//        virtual void neumVector(CoordSet&,ComplexVector&,
//                                double,double,double,double,int);
        virtual void neumVectorDeriv(CoordSet&,ComplexVector&,
                                     double,double,double,double,int,
                                     int pflag = 0);
        virtual void wetInterfaceVector(CoordSet&,ComplexVector&,
                                double,double,double,double,int,int);
        virtual void wetInterfaceVector(CoordSet&,ComplexVector&,
                                complex<double> (*)[3],complex<double>*);
        virtual void wetInterfaceVectorDeriv(CoordSet&,ComplexVector&,
                                complex<double> (*)[3],complex<double>*,
                                complex<double>*, int);
        virtual void sommerVector(CoordSet&, ComplexVector&, ComplexVector&);
        virtual void btVector(CoordSet&, ComplexVector&, ComplexVector&);

        virtual void ffpNeum(int, ComplexD*, CoordSet&,ComplexD*,
                             double,double(*)[3],double*);
        virtual void ffpDir(int, ComplexD*, CoordSet&,ComplexD*,ComplexD*,
                            double,double(*)[3],double*);
        virtual void ffpDir(int, ComplexD*, CoordSet&,ComplexD*,
                            double,double(*)[3],double*);
        virtual ComplexD ffpCoef(double k); 

	virtual void getNormal(CoordSet&, double[3]);
	virtual double getSize(CoordSet&);

        virtual FullSquareMatrixC turkelMatrix(CoordSet&, double, int);
        virtual FullSquareMatrixC turkelMatrix(CoordSet&, double, int, DComplex *);
        virtual bool isSommerElement() { return true; }
        virtual bool isPhantomElement() { return false; }

        virtual FullSquareMatrix stiffness(CoordSet&, double *d, int flg = 1);
        virtual FullSquareMatrix massMatrix(CoordSet&, double *mel, int cmflg=1);
        virtual FullSquareMatrix dampingMatrix(CoordSet&, double *cel, int cmflg=1);
        //virtual FullSquareMatrix imStiffness(CoordSet&, double *d, int flg = 1);

        virtual int ludcmp(FullSquareMatrix A, int n, int *indx);// LU decomposition
        virtual void lubksb (FullSquareMatrix A, int n, int *indx,double *b);//LU factorisation
        virtual void invert (FullSquareMatrix A, FullSquareMatrix B);

};

#endif

