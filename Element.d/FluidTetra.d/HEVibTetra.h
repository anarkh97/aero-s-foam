#ifndef _HEVIBTETRA_H_
#define _HEVIBTETRA_H_

#include <Element.d/Element.h>

class HEVibTetra: public Element {

	int nn[4];
public:
	HEVibTetra(int*);

	Element *clone();

	void renum(int *);

        FullSquareMatrix stiffness(CoordSet&, double *kel, int flg=1);
        FullSquareMatrix massMatrix(CoordSet&,double *mel, int cmflg=1);
	double           getMass(CoordSet& cs);
	void             markDofs(DofSetArray &);
        int*             dofs(DofSetArray &, int *p=0);
        int              numDofs();

        int              numNodes();
        int*             nodes(int * = 0);
        int 		 getTopNumber();

        bool             isHEVFluidElement()  { return true; }

	PrioInfo examine(int sub, MultiFront *)   {
           fprintf(stderr,"HEVibTetra.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C\n"); return *(new PrioInfo);
        };
        //void getThermalForce(CoordSet &, Vector &, Vector &force, int, GeomState *geomState=0) { force.zero(); }
private:

	double J[3][3];
	double Jinv[3][3];
	double dOmega;
	double gN[4][3];


	//double		volume(CoordSet&);
	double		volume() {return dOmega/6.0;} 
	void		computeMetrics(CoordSet&);
	void 		buildTetraStiff(double m[4][4]);
	
};
#endif

