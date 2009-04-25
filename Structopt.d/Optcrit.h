#ifndef _OPTCRIT_H_
#define _OPTCRIT_H_

#ifdef STRUCTOPT

#include <stdio.h>
#include <stdlib.h>

class Structopt;
struct critdata;
struct dldata;

//------------------------------------------------------------------------------
//                              Optimization Criteria 
//                          Base Class and Derived Classes
//                                                  created  9/1/98 by Kurt
//------------------------------------------------------------------------------

class Criteriatyp {
       
     public:

       //  Criteria-Type : 0 - Strain Energy
       //                  1 - Structural Mass
       //                  2 - Frequency
       //	           3 - Nodal Stress
       //	           4 - Nodal Displacement   
       //	           5 - Nodal Velocition   
       //	           6 - Nodal Acceleration   
       //	           7 - Kinetic Energy  
       //	           8 - Dissipated Energy   
       //	           9 - Aero Force / Sonic Boom  
       //	          10 - Elastic Force   
       //	          11 - Control Cost   
       //                 12 - StressIntegral
       //                 13 - DisplacementLeveling
       //                 14 - FORM-L2-norm
       //                 15 - Failure Probability
       //                 16 - Electrostatic
       //                 17 - Reliability Index
       //                 18 - Performance Measure
       //                 19 - Moment of Inertia
       //                 20 - Stress Leveling
       //                 21 - Nodal Temperature (from thermal)
       //                 22 - EM eigenvalue (pull-in)
       //                 23 - Change of Nodal Position
       //                 24 - Voltage (from electrostatic)
       //                 25 - Electric RHS (from electrostatic), ex. current
       //                 26 - Electrical Power
       //                 27 - Thermal Power Loss from Conduction
       //                 28 - Single Domain Temperature (non multi-physics)


       int typ;
       
     public:

        char * gettypname();

};

//------------------------------------------------------------------------------

class Stresstyp {
       
     public:

	//  Stress-Type : 00 -    von Mises on Lower  Surf.
	//                01 -    von Mises on Middle Surf.
	//                02 -    von Mises on Upper  Surf.
	//                10 - 1. Principal on Lower  Surf.
	//                11 - 1. Principal on Middle Surf.
	//                12 - 1. Principal on Upper  Surf.
	//                20 - 2. Principal on Lower  Surf.
	//                21 - 2. Principal on Middle Surf.
        //                22 - 2. Principal on Upper  Surf.

        int typ;

     public:

        char * gettypname();
};

//------------------------------------------------------------------------------

class Disptyp {
       
     public:

	//  Displacement-Type : 0 - x - Direction
	//                      1 - y - Direction
	//                      2 - z - Direction
	//                      3 - x - Rotation
	//                      4 - y - Rotation
	//                      5 - z - Rotation
        //                      6 - t - Length of Vektor

        int typ;

     public:

        char * gettypname();
};

//------------------------------------------------------------------------------

class Optcrit {

     public:	

        enum {StressAnalysis,ModalAnalysis,GeometricAnalysis,Explicit,
	      ExtCritAnalysis};

        int num;	    
        int anaId;
 
        double time;

        double val;
	double grad;

        Criteriatyp crittyp;           //Type of Criteria

     public:
	virtual ~Optcrit() {}
        static Optcrit * build (int, critdata&, double&, int);
	     
        int getnumber    () { return num ; }
        int gettyp       () { return crittyp.typ ; }	

        int checkTime(double,double); 
        int getAnalysisType();

        int getAnalysisId() { return anaId; };

        double getTime () { return time; }
	
	double getval  () { return val ;  }
	double getgrad () { return grad ; }
	
	void putval  ( double _val ) { val  = _val  ; }
	void putgrad ( double _grad) { grad = _grad ; }
        void setTime ( double _time) { time = _time ; }

        virtual int  getloc1()=0;
	virtual int  getloc2()=0;

        virtual double getRef() { return 0; }
        virtual double getPow() { return 0; }

        virtual int getEigvalId();
	virtual int getAeroforce (double &aforce) { return -1 ;}

        virtual void evaluate(Structopt*, double tmax=-1, double tmin=-1)=0;

        virtual void gradcrit(Structopt*, int)=0;
        virtual void gradpart(Structopt*)=0;
        virtual void graddu  (Structopt*)=0;
	
        virtual void print(FILE*)=0;
};

//------------------------------------------------------------------------------

class OptcritGlobal : public Optcrit {
    
     public:
     
        int    iloc1;
        int    iloc2;

        int    qFlag;
        int    nodeleTyp;

        int*   list;
        int    listsize;
 
        double refVal;
        double powFac;

        Stresstyp strtyp;
        Disptyp   distyp;	

     public:
     
        // strain energy, mass, frequency, aeroforces, control costs

        OptcritGlobal ( int _iloc1=0, int*  _list=0, int _listsize=0 );

        // mass in deformed configuration

        OptcritGlobal ( double& _d0, int* _list, int _listsize);

        // stress integral

        OptcritGlobal ( int _strtyp, double _d1, double _d2, int _i1, int _i2, 
                        int* _list, int _listsize); 

        int getloc1() { return iloc1; }
	int getloc2() { return iloc2; }

        double getRef() { return refVal; }
        double getPow() { return powFac; }

        int getEigvalId();
 	int getAeroforce (double &aforce);

	void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1);
	void gradcrit(Structopt *structopt,int);
	void gradpart(Structopt *structopt);
        void graddu  (Structopt *structopt);
	
        void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritNodstr : public Optcrit {
	
     public:
	        
	int    node;                                      //Node Number
        Stresstyp strtyp;                                 //Stress-Type

     public:
		     
        OptcritNodstr ( int _node = 0, int _strtyp = 0)
                      { node=_node; strtyp.typ=_strtyp; } 

        int getloc1() { return node; }
	int getloc2() { return strtyp.typ; }
                
        void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1);
        void gradcrit(Structopt *structopt,int);
        void gradpart(Structopt *structopt);
        void graddu  (Structopt *structopt);
	void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritDisp : public Optcrit {
	
     public:

	int node;                                       //Node Number

        Disptyp disptyp;                                //Displacement-Type
     
     public:
		     
        OptcritDisp ( int _node = 0, int _disptyp = 0)
                    { node        = _node; 
                      disptyp.typ = _disptyp; } 

        int getloc1() { return node; }
	int getloc2() { return disptyp.typ; }

        void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1); 
        void gradcrit(Structopt *structopt,int); 
        void gradpart(Structopt *structopt); 
        void graddu  (Structopt *structopt);
        void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritDispLevel : public Optcrit {
	
     public:

	int aFlag;               // averaging flag
        int qFlag;               // flag whether or not dividing by number nodes            

        int size;                // number of nodal displacements to be leveled
    
        int * nodeid;            // node id-number 
        int * distyp;            // type of displacement

        double * refVal;         // reference Value in nominator
        double * difVal;         // reference Value in denominator

        double   powFac;         // exponent
 
     public:
		     
        OptcritDispLevel ( int*, double *, dldata *, int);

        int getloc1() { return 0; }
	int getloc2() { return 0; }

        void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1); 
        void gradcrit(Structopt *structopt,int); 
        void gradpart(Structopt *structopt); 
        void graddu  (Structopt *structopt);
        void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritStressLevel : public Optcrit {
	
     public:

	int aFlag;               // averaging flag
        int qFlag;               // flag whether or not dividing by number nodes            
        int eFlag;               // flag whether elemental stress (1) or nodal stress (0)            

        int size;                // number of stresses to be leveled
    
        int * nodeid;            // node id-number 
        int * strtyp;            // type of stress

        double * refVal;         // reference Value in nominator
        double * difVal;         // reference Value in denominator

        double   powFac;         // exponent
 
     public:
		     
        OptcritStressLevel ( int*, double *, dldata *, int);

        int getloc1() { return 0; }
	int getloc2() { return 0; }

        void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1); 
        void gradcrit(Structopt *structopt,int); 
        void gradpart(Structopt *structopt); 
        void graddu  (Structopt *structopt);
        void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritExt : public Optcrit {
 
    public:

	int node;                                       //Node Number
	int getloc2() { return 0; }
	
        OptcritExt( int _node = 0 ){ node=_node;} 
        int getloc1() { return node; }
        void evaluate(Structopt *structopt, double tmax=-1, double tmin=-1); 
        void gradcrit(Structopt *structopt, int); 
        void gradpart(Structopt *structopt); 
        void graddu  (Structopt *structopt);
        void print(FILE*);
};

//------------------------------------------------------------------------------

class OptcritElecMech : public Optcrit {
 
    public:

        OptcritElecMech() {} 

        int  getloc1() { return 0; }
	int  getloc2() { return 0; }

        void evaluate(Structopt *structopt,double tmax=-1,double tmin=-1); 
        void gradcrit(Structopt *structopt,int); 
        void gradpart(Structopt *structopt); 
        void graddu  (Structopt *structopt);
        void print(FILE*);
};
#endif

#endif
