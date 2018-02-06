#ifndef _DEMHELM2D_H_
#define _DEMHELM2D_H_

#include <cmath>
#include <Element.d/DEM.d/DEMElement.h>

class DGMHelm2d_LM: public DEMLM {
public:
	void init() override {};
	virtual complex<double> coef(int)=0;
};


class DGMHelm2d_1_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 1; }
	int type() const override { return 1; }
	complex<double> coef(int) override;
};

class DGMHelm2d_2_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 2; }
	int type() const override { return 2; }
	complex<double> coef(int) override;
};

class DGMHelm2d_4_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 4; }
	int type() const override { return 3; }
	complex<double> coef(int) override;
};

class DGMHelm2d_8_LM: public DGMHelm2d_LM {
public:
	int nDofs() const override { return 8; }
	int type() const override { return 4; }
	complex<double> coef(int) override;
};

class DGMHelm2d_Eva_LM: public DEMLM {
public:
	virtual complex<double> ldir(int,double[2])=0;
};

class DGMHelm2d_Eva1_2_LM: public DGMHelm2d_Eva_LM {
	int ndofs;
	double gamma;
	double normal[2];
public:
	DGMHelm2d_Eva1_2_LM() {
		// needs to be initialized
		ndofs = 0;
	}
	virtual void init();
	virtual int nDofs() const { return ndofs; }
	int type() const override { return 5; }
	virtual complex<double> ldir(int,double[2]);
};


class DGMHelm2d: public DEMElement {
public:
	int ndir;
	virtual complex<double> dir(int);
	int o;
	DGMHelm2d(int _o, int* nodenums);
	int defaultLMType() const override { return 1; }
	bool dgmFlag() const override { return true; }
	bool condensedFlag() const override  {
		for(int i=0;i<nFaces();i++) if (bc[i]==3) return false;
// return false;
		return true;
	}
// bool condensedFlag() const override  { return false; }
	virtual bool storeMatricesFlag() { return true; }
// virtual bool storeMatricesFlag() { return false; }
	int nPolynomialDofs() const  { return 0; }
	virtual int nEnrichmentDofs() const override { return ndir; }
	int nGeomNodes() const override { return (o>0)?o*o:((-o)*(-o+1))/2; }
	int nFaces() const override { return (o>0)?4:3; }
	int nFaceCorners(int fi) const override { return 2; }
	virtual int *faceCorners(int fi) {
		int *fc = new int[2];
		if (o>0) {
			if (fi==1) { fc[0] = nn[0]; fc[1] = nn[o-1]; }
			else if (fi==2) { fc[0] = nn[o-1]; fc[1] = nn[o*o-1]; }
			else if (fi==3) { fc[0] = nn[o*o-1]; fc[1] = nn[o*(o-1)]; }
			else if (fi==4) { fc[0] = nn[o*(o-1)]; fc[1] = nn[0]; }
			return fc;
		} else {
			if (fi==3) { fc[0] = nn[0]; fc[1] = nn[-o-1]; }
			else if (fi==1) { fc[0] = nn[-o-1]; fc[1] = nn[((-o)*(-o+1))/2-1]; }
			else if (fi==2) { fc[0] = nn[((-o)*(-o+1))/2-1]; fc[1] = nn[0]; }
			return fc;
		}
	}
	virtual int polyDofType() const override { return DofSet::Helm; }
	virtual int polyDofsPerNode() const override { return 1; }

	virtual void getRef(double *xyz,double *xy);
	virtual void createM(complex<double>*);
	virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
	virtual void createRHS(complex<double>*);
	void createSol(double *xyz, complex<double>*,
	               complex<double>*) override;
};


class DGMHelm2d_4: public DGMHelm2d {
public:
	DGMHelm2d_4(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 1; }
};

class DGMHelm2d_4t: public DGMHelm2d {
public:
	DGMHelm2d_4t(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 1; }
};

class DGMHelm2d_8: public DGMHelm2d {
public:
	DGMHelm2d_8(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 2; }
};

class DGMHelm2d_8t: public DGMHelm2d {
public:
	DGMHelm2d_8t(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 2; }
};

class DGMHelm2d_16: public DGMHelm2d {
public:
	DGMHelm2d_16(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 3; }
};

class DGMHelm2d_32: public DGMHelm2d {
public:
	DGMHelm2d_32(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 4; }
};


class DGMHelm2d_Eva: public DGMHelm2d {
public:
	double normal[2];
	double gamma;
	DGMHelm2d_Eva(int _o, int* nodenums);
};

class DGMHelm2d_Eva2_8: public DGMHelm2d_Eva {
public:
	DGMHelm2d_Eva2_8(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 5; }
};


class DEMHelm2d: public DGMHelm2d {
public:
	DEMHelm2d(int _o, int* nodenums);
	bool dgmFlag() const override { return false; }
	int nPolynomialDofs() const  { return (o>0)?o*o:((-o)*(-o+1))/2; }

	virtual void createM(complex<double>*);
	virtual void interfMatrix(int fi, DEMElement*,complex<double>*);
	virtual void createRHS(complex<double>*);
};


class DEMHelm2d_4: public DEMHelm2d {
public:
	DEMHelm2d_4(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 1; }
};

class DEMHelm2d_4t: public DEMHelm2d {
public:
	DEMHelm2d_4t(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 1; }
};

class DEMHelm2d_8: public DEMHelm2d {
public:
	DEMHelm2d_8(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 2; }
};

class DEMHelm2d_8t: public DEMHelm2d {
public:
	DEMHelm2d_8t(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 2; }
};

class DEMHelm2d_16: public DEMHelm2d {
public:
	DEMHelm2d_16(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 3; }
};

class DEMHelm2d_32: public DEMHelm2d {
public:
	DEMHelm2d_32(int _o, int* nodenums);
	virtual complex<double> dir(int);
	int defaultLMType() const override { return 4; }
};

#endif

