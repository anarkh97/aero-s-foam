#ifndef	_ELEMENT_H_
#define	_ELEMENT_H_

#include <Math.d/matrix.h>
#include <Utils.d/BlockAlloc.h>
#include <Utils.d/dofset.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <complex>

// this is a fix to get around apparent template bug in solaris compiler
inline double abs(std::complex<double> a)
{
  return sqrt(a.real()*a.real() + a.imag()*a.imag());
}

class Corotator;
class State;
class PolygonSet;
class InterpPoint;
class Material;
class LMPCons;
class GeomState;  // PJSA: for sgi intel
template <class Scalar> class GenVector; // PJSA: for sgi intel
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;  // PJSA: for sgi intel
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
//template <class T> class ResizeArray;

// Boundary Condition Structure
struct BCond {
  int nnum;   // node number
  int dofnum; // dof number (0-6)
  double val;    // value of bc
  void setData(int _nnum, int _dofnum, double _val) { nnum = _nnum; dofnum = _dofnum; val = _val; };
};

// Complex Boundary Condition Structure
struct ComplexBCond {
  int    nnum;   // node number
  int    dofnum; // dof number
  double reval;  // real value of bc
  double imval;  // imaginary value of bc
  void setData(int _nnum, int _dofnum, double _reval, double _imval) { nnum = _nnum; dofnum = _dofnum; reval = _reval; imval = _imval; };
};


struct NumExc {
  int gN, n1, n2, n3;
  NumExc(int g, int a1, int a2, int a3)
     { gN = g; n1 = a1; n2 = a2; n3 = a3; }
};

typedef double EFrame[3][3];

// ****************************************************************
// * Class StructProp contains the defined structural properties  *
// ****************************************************************

// contains material and geometrical properties of elements

struct PMLProps {
 int PMLtype;
 double gamma, Rx, Ry, Rz, Sx, Sy, Sz;
};

class StructProp {
  public:
    union
        {
	double	E;      // Elastic modulus
	double d0;      // Initial stiffness of nonlin element
        double ky;
	};
    union {
	double  A;  	// Cross-sectional area
        double kx;
        };
    union
	{
	double	nu; 	// Poisson's ratio
	double a;	// shear-t contribution ratio
	double Ep;      // Plastic modulus
        double kz;
        double lambda;  // damage control
	};

	double  rho; 	// Mass density per unit volume
     union
	{
	double  eh;	// Element thickness
	double  C1;	// Non-uniform torsion constant
	double b;       // shear-s contribution ratio
	double xo;      // Plastic parameter    -January 2002 - JMP
	};
     union
        {
	double  Ixx;	// Cross-sectional moment of inertia about local x-axis
        double  ss;     // speed of sound
        };
	double  Iyy;	// Cross-sectional moment of inertia about local y-axis
	double  Izz;	// Cross-sectional moment of inertia about local z-axis
     union {
	double c; 	// Thermal convection coefficient
	double alphaY;  // Shear deflection constant associated to Iyy
	double sigmax;
	double sigE;    // Elastic limit
	double ft;      // Tensile strength     -October 2001 - JMP
        double eps;     // Radiation: Emissivity of the body (for a black body, eps=1)
	};
     union {
	double k;       // Heat conduction coefficient
	double alphaZ;  // Shear deflection constant associated to Izz
	double v2;      // Fracture Energy
	double fc;      // Compressive strength -October 2001 - JMP
	};
     union {
        double Q;	// Specific heat coeffiecient
        double betaDamp; // Rayleigh stiffness damping coefficient
        };
     union {
	double W;	// Thermal expansion coefficient
        double alphaDamp; // Rayleigh mass damping coefficient
        };
        double P;	// Perimeter/circumference of thermal truss/beams
     union {
        double Ta;	// Ambient temperature
        double Tr;      // Temperature of the enclosure receiving the radiation
        };
        double sigma;   // Stefan's constant (5.6704e-8 in SI)
	double ymin;    // minimum height (< 0) for cross section of beam (local y-direction)
	double ymax;    // maximum height (> 0) for cross section of beam (local y-direction)
	double zmin;    // minimum height (< 0) for cross section of beam (local z-direction)
	double zmax;	// maximum height (> 0) for cross section of beam (local z-direction)

        double kappaHelm; // wave number for Helmholtz proplem
        double kappaHelmImag; // imaginary part of the wavenumber for
                              // Helmholtz problem
        complex<double> soundSpeed;

	// Fabric Material Options
	int F_op; // Fabric Material Option
	double F_Uc; // Critical Stretch of the Fibrils
	double F_Uf; // Failure Stretch of the Yarn
     union {
	double F_h; // Initial Height of the Yarn
	double F_Estd; // Standard Deviation in E
	};
     union {
	double F_d; // Standard Deviation for the Fibril Inclination
	double F_lamstd; // Standard Deviation in Lambda
	};
	double F_dlambda; // Change in Lambda as a function of E
	int F_np; // Number of Points in Curve Fit for Lambda
	int F_Nf; // Number of Fibrils in a Yarn
	int Seed; // Seed for Random Number Generator


        PMLProps fp;

        bool isReal;

	/** the W and E coefficient might encode integer values when they're negative
	 * (see manual for this). Heavily templated Sower needs a temporary storage that's adressable.
	 * For some reason it wont let us do that any other way than this. I'm hoping to improve this. TG
	 */
	int __SOWER_TMP;

        StructProp() { E = 0.0; A = 0.0; nu = 0.0; rho = 1.0; eh = 0.0; Ixx = 0.0;
                       Iyy = 0.0; Izz = 0.0; c = 0.0; k = 0.0; Q = 0.0; W = 0.0;
                       P = 0.0; Ta = 0.0; sigma = 0.0;
                       kappaHelm = 0.0; kappaHelmImag = 0.0; fp.PMLtype = 0;
                       soundSpeed = 1.0;
                       ymin = 0.0; ymax = 0.0;
		       zmin = 0.0; zmax = 0.0; isReal = false; } // Undefined property

        StructProp(StructProp *p)  {E=p->E; A=p->A; nu=p->nu; rho=p->rho; eh=p->eh;
                                    Ixx=p->Ixx; Iyy=p->Iyy; Izz=p->Izz; c=p->c;
                                    k=p->k; Q=p->Q; W=p->W; P=p->P; Ta=p->Ta; sigma=p->sigma;
                                    ymin=p->ymin; ymax=p->ymax; zmin=p->zmin; zmax=p->zmax;
                                    kappaHelm = p->kappaHelm;
                                    kappaHelmImag = p->kappaHelmImag;
                                    soundSpeed = p->soundSpeed;
                                    fp.PMLtype = p->fp.PMLtype; isReal = p->isReal; }

        //StructProp(const StructProp &p) { cerr << "in StructProp copy constructor\n"; }
        //StructProp& operator = (const StructProp &p) { cerr << "in StructProp  assignment operator\n"; }
};

// ****************************************************************
// *
// *     class Node: Keeps the coordinates of a node
// *
// ****************************************************************

class Node {
  public:
	// Constructors
        Node() {}
        Node(double *xyz) { x = xyz[0]; y = xyz[1]; z = xyz[2]; }
        Node(const Node &node) { x = node.x; y = node.y; z = node.z; }
	~Node() {};
	double distance2(const Node& node) const
	   { return (node.x-x)*(node.x-x)+(node.y-y)*(node.y-y)+(node.z-z)*(node.z-z); }
        double distance(const Node &node) const { return sqrt(distance2(node)); }

	// Coordinates
	double 	x;
	double 	y;
	double 	z;
};

// ****************************************************************
// *                        WARNING                               *
// *       Nodes in the CoordSet are from 0 and not from 1        *
// *       It is a good idea to subtract 1 from node input        *
// *       from the user at any time a node is read in rather     *
// *       then keeping it from 1 and then subtracting 1 here     *
// *       and there....                                          *
// *       functions for this class are in Element.d/Elemset.C   *
// ****************************************************************

class CoordSet {
        int nmax;
	Node **nodes;
        BlockAlloc ba;
public:
        // Constructors
        CoordSet(int = 256);
	CoordSet(const CoordSet&, int n, int *nd);
        ~CoordSet() { deleteNodes(); }

	CoordSet* copy();

	// Member functions
        int size() { return nmax ; }
        int last();
	void list(void)
	  {
	    std::cerr << "Size: " << nmax << std::endl;
	    for(int i = 0;i<nmax-1;i++)
	      if(nodes[i])
		std::cout << i << " : " << nodes[i]->x << "," << nodes[i]->y << "," << nodes[i]->z << std::endl;
	      else
		std::cerr << "warning no assignation for node #" << i <<std::endl;
	  }
        int exist(int numNodes);
        void  nodeadd(int n, double*xyz);
        void  nodeadd(int n, Node &node);
        void  nodeCopy(int n, CoordSet &nd);
	Node &getNode(int n);
        void getCoordinates(int *nn, int numNodes,
                            double *xx, double *yy, double *zz);

        Node *operator[] (int i) { return (i >= nmax) ? 0 : nodes[i]; }
        Node *operator[] (int i) const { return (i >= nmax) ? 0 : nodes[i]; }

        void deleteNodes()
	  { if(nodes) {
	      delete [] nodes; nodes = 0; nmax = 0;
	    }
	  }
};


 //DEC
struct PrioInfo {
  int priority; // the level of priority (undefined in isReady == 0)
  bool isReady; // whether this element is ready to be inserted
  // PrioInfo() { isReady = false; }  PJSA
};

class MultiFront;
// END DEC


/** \brief Abstract Element class
 *
 * Class Element defines functions for finite elements.         *
 * Each element has a structural property and a pressure        *
 * associated with it. Each element defines it's own dofs, sets *
 * it's own node numbers, calculates stiffness, calculates mass *
 * von mises stress, internal force and displacements. Where    *
 * these functions are written for each type of element and if  *
 * a function is not defined for an element type, it will       *
 * default to the appropriate function in Element.C             *
 * i.e. an element with a zero mass matrix, will just return a  *
 * matrix of the appropriate size containing zeroes.            *
 ****************************************************************/

class Element {
  public:
        enum Category { Structural, Acoustic, Thermal, Sloshing, HEV, Undefined };
  private:
        Category category;
        double _weight, _trueWeight;
        int elementType;
  protected:
	StructProp *prop;	// structural properties for this element
        double pressure;	// pressure force applied to element
        double preload;
        bool myProp;
        int glNum;
        vector<double> factors;
	void lumpMatrix(FullSquareMatrix&);
  public:
        Element() { prop = 0; pressure = 0.0; preload = 0.0;
        _weight=1.0; _trueWeight=1.0; myProp = false; category = Undefined; };
        virtual ~Element() { if(myProp && prop) delete prop; }
        StructProp * getProperty() { return prop; }

        virtual Element *clone() { return 0; }
        virtual void renum(int *)=0;

        static Element *build(int,int,int*);

        virtual void buildCorotator(CoordSet &)  {};

        virtual void setProp(StructProp *p, bool _myProp = false) {
          if(myProp && prop)  {
            delete prop;
            prop=0;
          }
          prop = p; myProp = _myProp;
        }
	virtual void setPressure(double pres) { pressure = pres; }
        virtual double getPressure() { return pressure; }

        // By default ingore any element preload
        virtual void setPreLoad(double load, int &flg) { }
        virtual double getPreLoad() { return preload; }

        void setGlNum(int gn) { glNum = gn; }
        int getGlNum()  { return glNum; }
        virtual void setFrame(EFrame *) {} // By default ignore the frame
        // By default an element does not need a frame
        virtual int buildFrame(CoordSet&) { return 1; }
        virtual void setOffset(double*){}
	    virtual void setCompositeData(int _type, int nlays, double *lData,
                                      double *coefs, double *frame);
        virtual double * setCompositeData2(int _type, int nlays, double *lData,
                                           double *coefs, CoordSet &cs, double theta);

        virtual FullSquareMatrix stiffness(CoordSet& cs,double *k,int flg=1);
	virtual FullSquareMatrix massMatrix(CoordSet& cs,double *m,int cmflg=1);
        virtual FullSquareMatrix imStiffness(CoordSet& cs,double *k,int flg=1);
	FullSquareMatrix massMatrix(CoordSet& cs, double* m, double mratio);
        virtual FullSquareMatrixC stiffness(CoordSet&, complex<double> *d) {return FullSquareMatrixC();};    //CRW
        virtual FullSquareMatrixC massMatrix(CoordSet&, complex<double> *d) {return FullSquareMatrixC();};    //CRW

	virtual FullSquareMatrix dampingMatrix(CoordSet& cs,double *m,int cmflg=1);

        virtual double getMass(CoordSet&) { return 0; }
        virtual double getDCmass(CoordSet &,Vector &, double&) { return 0; }

        virtual void   getGravityForce(CoordSet&,double *gravity,Vector &force,
                                       int gravflg, GeomState *gs=0);

        virtual void   getThermalForce(CoordSet& cs,Vector &ndT,Vector &force,
                                       int glflag, GeomState *gs=0);

	virtual void   getIntrnForce(Vector &elForce, CoordSet& cs,
				     double *elDisp, int Index, double *ndTemps);

        // this can't be templated, c++ doesn't allow virtual member functions to be templated
	virtual void   getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
		                   Vector &elDisp, int strInd, int surface=0,
                                   double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum = 0); //CBM

        virtual void   getVonMises(ComplexVector &stress, Vector &weight, CoordSet &cs,
                                   ComplexVector &elDisp, int strInd, int surface=0,
                                   double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum = 0); //CBM

	virtual void   getVonMisesInt(CoordSet &,Vector &,double &,double &, int,
				      double &,double &, double* dT=0 );

        virtual void   getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                                    Vector &elDisp, int strInd, int surface=0,
                                    double *ndTemps=0);

        virtual void   computeHeatFluxes(Vector& heatflux, CoordSet &cs,
                                         Vector& elTemp, int hflInd);

        virtual void  trussHeatFluxes(double &trussflux, CoordSet &cs,
                                      Vector& elTemp, int hflInd) {}

        //ADDED FOR SLOSHING PROBLEM, EC, 20070723
        virtual void   computeSloshDisp(Vector& fluidDispSlosh, CoordSet &cs,
                                         Vector& elPotSlosh, int hflInd);

        //ADDED FOR SLOSHING PROBLEM, EC, 20081101
        virtual void   computeSloshDispAll(Vector& fluidDispSlosh, CoordSet &cs,
                                         Vector& elPotSlosh);

        virtual void   computeDisp(CoordSet&, State&, const InterpPoint &,
                                   double*, GeomState *gs=0)  {}

        virtual void   getFlLoad(CoordSet &, const InterpPoint &,
                                 double *, double *, GeomState *gs=0) {}

        virtual void   computeTemp(CoordSet&, State &, double[2],
                                   double*)  {}
        virtual void   getFlFlux(double[2], double *, double *) {}

	virtual void   markDofs(DofSetArray &)= 0;
	virtual int*   dofs(DofSetArray &, int *p=0)=0;
	virtual int    numDofs()=0;

	virtual int    numNodes() = 0;
        virtual int    totNumNodes() { return numNodes(); }
	virtual int*   nodes(int * = 0) = 0;
        virtual int*   allNodes(int *x = 0) { return nodes(x); }

        virtual Corotator *getCorotator(CoordSet &, double *, int = 2, int = 2)
        { printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0; }

	virtual int  getTopNumber();
        virtual int  numTopNodes() { return numNodes(); }   // this is the number of nodes printed in the top file
                                                            // can make it different to numNodes for elements that aren't
                                                            // supported by xpost eg RigidSolid6Dof
        virtual void computePressureForce(CoordSet& cs,Vector& elPressureForce,
                                          GeomState *gs=0, int cflg = 0);
        virtual void computePressureForce(Node *nd,Vector& elPressureForce,
                                          double *gs=0, int cflg = 0);
	virtual double * getMidPoint(CoordSet &)  { return 0; }
	/* toxa: use this midPoint instead */
	virtual void getMidPoint(CoordSet &, double* result)  { assert(0); }
        virtual double * getCompositeData(int )   { return 0; }
        virtual double * getCompositeFrame()      { return 0; }

        virtual int      getCompositeLayer()      { return 0; }

	virtual double   getMoment(Vector&, CoordSet&, int, int)
	                 { return 0.0; }

        virtual int     dim() { return 3; }

        virtual void addFaces(PolygonSet *pset);

        virtual void setMaterial(Material *);

	virtual int numInternalNodes() { return 0; }

	virtual void setInternalNodes(int *) {}

	virtual bool isSafe() { return true; }


	// from DEC
	// TOP/DOMDEC Element Functions
	virtual int facelist(PolygonSet &, int * = 0) {return 0; }


	// DEC : Routines for the decomposer
	// isStart indicates if an element is suitable to
	// be a germination center for a subdomain (bars are not)
	virtual bool isStart() { return true; }

	virtual bool isSpring() { return false; }

	virtual bool isMass() { return false; }

	virtual bool hasRot() { return false; }
 	virtual PrioInfo examine(int sub, MultiFront *mf);
	virtual double weight() { return _weight; }
	virtual double trueWeight() { return _trueWeight; }
	void setWeight(double weight){ _weight = weight; }
	void setTrueWeight(double trueWeight){ _trueWeight = trueWeight; }

	void getCG(CoordSet &cset, double &xcg, double &ycg, double &zcg);
	// END FROM DEC

         // PJSA: this need to be defined for 6 node tri shell & 8 node quad shell
        virtual bool isRotMidSideNode(int iNode) { return false; }

	virtual bool hasDamping() { return false; }
        bool isFluidElement();
        virtual bool isSommerElement() { return false; }
        virtual bool isMpcElement() { return false; }
        virtual bool isFsiElement() { return false; }
        virtual bool isHEVFluidElement() { return false; }  //ADDED FOR HEV PROBLEM, EC, 20070820
        virtual int  fsiFluidNode() { return -1; }
        virtual int  fsiStrutNode() { return -1; }
        virtual bool isRigidMpcElement(const DofSet & = DofSet::nullDofset, bool forAllNodes=false) { return false; }
        virtual bool isRigidElement() { return false; }
        virtual void computeMPCs(CoordSet &cs) { }
        virtual int getNumMPCs() { return 0; }
        virtual LMPCons** getMPCs() { return 0; }

        virtual double helmCoef() { return prop->kappaHelm * prop->kappaHelm; }
        virtual complex<double> helmCoefC() { return
           complex<double>(prop->kappaHelm,prop->kappaHelmImag)*
           complex<double>(prop->kappaHelm,prop->kappaHelmImag);
        }
        virtual double helmCoef(double omega) {
           return omega*omega*real(prop->soundSpeed)*real(prop->soundSpeed);
        }
        virtual complex<double> helmCoefC(double omega) {
           return (omega*omega)*prop->soundSpeed*prop->soundSpeed;
        }

        virtual bool isShell() { return false; }

        bool isConstraintElement() { return (isRigidElement() || isMpcElement() || isRigidMpcElement() || isFsiElement()); }
        virtual bool isPhantomElement() { return (!(prop || isConstraintElement() || isSommerElement())); }

        int getElementType() { return elementType; }
        void setElementType(int _type) { elementType = _type; }

	// friend class Domain;
        friend class SuperElement;

	// complex interface for the Element
	virtual bool isComplex() { return false; }
        virtual FullSquareMatrixC complexStiffness(CoordSet& cs, DComplex *k,int flg=1);
	virtual FullSquareMatrixC complexDampingMatrix(CoordSet& cs, DComplex* c, int flg=1);
	virtual FullSquareMatrixC complexMassMatrix(CoordSet& cs, DComplex* m, double mratio);

        Category getCategory() { return category; } 
        void setCategory(Category _category) { category = _category; } // currently this is only called for thermal elements, could be extended.
        bool isDamped() { return (getCategory() != Thermal && !isSpring()) ? (prop && (prop->alphaDamp != 0.0 || prop->betaDamp != 0.0)) : false; }

        virtual int getMassType() { return 1; }  // 0: lumped, 1: consistent, 2: both
                                                 // notes: (a) if getMassType returns 0 then lumped gravity force will always be used for dynamics
                                                 //        (b) is getMassType returns 1 then lumping is done using diagonal scaling if required (default)
};

// *****************************************************************
// *                        WARNING                                *
// *       The same remark as for node is valid for elements       *
// *       The functions for this class are in Element.d/ElemSet.C *
// *****************************************************************

class Elemset
{
  protected:
    Element **elem;
    int emax;
    int nPhantoms;
    BlockAlloc ba;
    bool myData;
  public:
    Elemset(int = 256);
    virtual ~Elemset() { deleteElems(); }
    Elemset(Elemset &globalset, int numlocal, int *localToGlobalMap);
    int size() const { return emax; }
    int last() const;
    int numPhantoms()  { return nPhantoms; }
    Element *operator[] (int i) const { return elem[i]; }
    void elemadd(int num, Element *);
    void elemadd(int num, int type, int nnodes, int *nodes);
    void mpcelemadd(int num, LMPCons *mpc);
    void fsielemadd(int num, LMPCons *fsi);
    void setEmax(int max)  { emax = max; }
    void setNumPhantoms(int n)  { nPhantoms = n; }
    void list();
    void deleteElems()  { if(elem) { if(myData) for(int i=0; i<last(); ++i) delete elem[i]; delete [] elem; elem = 0; emax = 0;} }
    //void deleteElems()  { if(elem) delete [] elem; elem = 0; emax = 0; }
    void remove(int num) { elem[num] = 0; }//DEC
    void setMyData(bool _myData) { myData = _myData; }
    bool hasDamping() { for(int i=0; i<last(); ++i) if (elem[i]->isDamped()) return true; return false; }
   void collapseRigid6(); // replaces all 6-DOF rigid elements that share a node by a single element
};

class EsetGeomAccessor {
  public:
    static int getNum(/*const*/ Elemset &eSet, int elNum)
     { return eSet[elNum] ? eSet[elNum]->totNumNodes() : 0; }
    static int getSize(const Elemset &eSet)
     { return eSet.last(); }
    static int *getData(/*const*/ Elemset &eSet, int elNum, int *nd)
     { return eSet[elNum] ? eSet[elNum]->allNodes(nd) : 0; }
};

class EsetDataAccessor {
 public:
  static int getNum(/*const*/ Elemset &eSet, int elNum)
    { return 1; }
    static int getSize(const Elemset &eSet)
     { return eSet.last(); }
    static int *getData(/*const*/ Elemset &eSet, int elNum, int *nd)
     {
       if(eSet[elNum])
	 {
	   if(nd) { nd[0] = elNum; return nd;}
	   else { std ::cerr << "(EE) error : unefficient use of EsetPressureAccessor" << std::endl; return 0;}
	 }
       else
	 {
	   return 0;
	 }
     }
};

#include <Driver.d/EFrameData.h>

class  EFrameDataAccessor{
 public:
    static int getNum(/*const*/ std::pair<int,ResizeArray<EFrameData>* > &, int )
     { return 1; }
    static int getSize(const std::pair<int,ResizeArray<EFrameData>* > &o)
     { return o.first; }
    static int *getData(/*const*/ std::pair<int,ResizeArray<EFrameData>* > &o, int i, int *nd)
     {
       if(i>=o.first)
	 cout << "EFrameDataAccessor corruption" << endl;
       if(nd) { nd[0] = ((*(o.second))[i]).elnum; return nd; }
       else return &(*o.second)[i].elnum;
     }
};

#include <Driver.d/DMassData.h>

class DimassAccessor{
 public:
  static int getNum(/*const*/ vector<DMassData*> &, int )
    { return 1; }
  static int getSize(const  vector<DMassData*> &o)
    { return o.size(); }
  static int *getData(/*const*/  vector<DMassData*> &o, int i, int *nd)
    {
      if(nd != 0)
	{
	  nd[0] = (*(o[i])).node; return nd;
	}
      else
	return &(*(o[i])).node;
    }
};

class BCDataAccessor {
  public:
    static int getNum(/*const*/ std::pair<int,BCond *> &, int )
     { return 1; }
    static int getSize(const std::pair<int,BCond *> &o)
     { return o.first; }
    static int *getData(/*const*/ std::pair<int,BCond *> &o, int i, int *nd)
     {
       if(i>=o.first)
	 cout << "BCDataAccessor" << endl;
       if(nd) { nd[0] = o.second[i].nnum; return nd; }
       else return &o.second[i].nnum;
     }
};

class ComplexBCDataAccessor {
  public:
    static int getNum(/*const*/ std::pair<int,ComplexBCond *> &, int )
     { return 1; }
    static int getSize(const std::pair<int,ComplexBCond *> &o)
     { return o.first; }
    static int *getData(/*const*/ std::pair<int,ComplexBCond *> &o, int i, int *nd)
     {
       if(i>=o.first)
         cout << "ComplexBCDataAccessor" << endl;
       if(nd) { nd[0] = o.second[i].nnum; return nd; }
       else return &o.second[i].nnum;
     }
};

class ElementFactory
{
 public:
  ElementFactory() {}
  virtual ~ElementFactory() {}
  virtual Element* elemadd(int num, int type, int nnodes, int *nodes,
			   BlockAlloc& ba);
};

#endif
