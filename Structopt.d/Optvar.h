#ifndef _OPTVAR_H_
#define _OPTVAR_H_

#ifdef STRUCTOPT

#include <cstdio>
#include <Utils.d/resize_array.h>
#include <vector>

class Structopt;
class Optopr;
class OptActInfo;

struct funcall;
struct absvardata;

//------------------------------------------------------------------------------

class Absvar {
 protected:
  int        num;             //Number of Abstract Variable
  
  int        typ;             //Type of Abstract Variable
  // 0 : design variable
  // 1 : random variable
  
  double     val;             //Value of AbstractVariable     
  double     grad;            //Gradient of AbstractVariable
  double     scal;            //Scaling Factors of Variable 
  double     upp;             //Upper Bounds on Variable
  double     low;             //Lower Bounds on Variable
  
  double     newval;          //Value of new time step in dynamics
  double     dynscl;          //Blending coefficient in dynamics
  
 public:
  virtual ~Absvar() {}
  static Absvar* create(int,absvardata&);
  
  int  getnumber()              { return num; }
  int  gettype()                { return typ; }
  
  double getval ()              { return val; }
  double getgrad()              { return grad; }
  double getlow()               { return low; }
  double getupp()               { return upp; }
  double getscal()              { return scal; }
  double getdynscl()            { return dynscl; }
  void   setdynscl(double newscal) { dynscl = newscal; }
  double getnewval()            { return newval; }
  
  void setnewval()              { double dum = newval; newval=val; val=dum; }
  void saveval()                { newval = val; }
  void setval(double nval)      { val    = nval; }
  
  
  // virtual functions required for all derived classes
  virtual void   print(FILE*)=0;

  virtual double getScaledValue()=0;     
  virtual double getScaledGrad()=0;     
  virtual double getScaledUpperBound()=0;  
  virtual double getScaledLowerBound()=0;  
  
  virtual void putScaledValue(double&)=0; 
  virtual void setgrad(double&)=0;        
  
  // functions related only to random variables  
  virtual int    getDistribution();
  virtual double getMeanValue(); 
  virtual double getStdDev();
  virtual void   generateRandom();
  virtual void   setMeanValue(double);
  virtual void   resetMeanValue();
  virtual void   setFormFlag(int);
  virtual int    getFormFlag();
};

//------------------------------------------------------------------------------


class Dsgvar:public Absvar {

   public:
   
     Dsgvar(absvardata&);

     void   print(FILE*);

     double getScaledValue()        { return val*scal;  }
     double getScaledGrad()         { return grad*scal; }
     double getScaledUpperBound()   { return upp*scal;  }
     double getScaledLowerBound()   { return low*scal;  }

     void putScaledValue(double& x) { val = x/scal;  }
     void setgrad(double& g)        { grad = g/scal; }

};

//------------------------------------------------------------------------------

class Rndvar:public Absvar {

   public:

     enum {normal,uniform,lognormal};
   
     int    dist;	     //Type of statistical distribution
     int    formFlag;	     //Flag indicating FORM transformation required

     double mean;	     //Mean Value of Random Variable	 
     double sdev;	     //Standard Deviation     
   
   public:
   
     Rndvar(absvardata &);

     void print(FILE*); 

     double getScaledValue();
     double getScaledGrad();
     double getScaledUpperBound();
     double getScaledLowerBound();

     void putScaledValue(double&);
     void setgrad(double&);

     int    getDistribution()      { return dist; }
     double getMeanValue()         { return mean; } 
     double getStdDev()            { return sdev; }
     void   setMeanValue(double v) { val=v;    }
     void   resetMeanValue()       { val=mean;    }
     void   generateRandom();
     void   setFormFlag(int);
     int    getFormFlag()     { return formFlag; }

};

//------------------------------------------------------------------------------

class Stcvar {
 protected:
     int rflag;                     //Indicates use in design/reliability

     int num;                       //Number of Structural Variable   
     int typ;                       //Type of Structural Variable
     
     int loc1;                      //Specification
     int loc2;
     
     double val;                    //Value
     double orgval;                 //Original Value                
     double grad;                   //Gradient
     
     double *pval;                  //Pointer to structural Quantity
     double *pgrad;                 //Pointer to Gradient

   public:
     Stcvar() : rflag(0), num(0), typ(0), loc1(0), loc2(0), pval(0), pgrad(0) {}
     static Stcvar * build( int _num=0,  int _typ=0, 
                            int _loc1=0, int _loc2=0, int rfg=0); 

     int  getnumber()             { return num ; }
     int  gettype() { return typ; }

     double getval() { return val; }
     double getgrad() { return grad; }
     
     int  getDesignVelocityComp() { int rval=-1; 
                                    if (typ == 2) rval=loc2; 
                                    return rval; }
				
     int getGroupType()           { return typ ; }			

     char * gettypname();
    
     void setval  ( double _val )  { val  = _val ; }
     void setgrad ( double _grad ) { grad = _grad ; }
    
     void updval();
     void updgrad();                
     void resetval();

     void setvalptr  (Structopt *structopt);
     void setgradptr (Structopt *structopt);

     void print(FILE*);

};

//------------------------------------------------------------------------------

class Optvar {
 private:
  std::vector<Absvar*>   absvar;      //Abstract Variables - Random Variables 
  std::vector<Stcvar*>   stcvar;      //Structural Variables     
  std::vector<Optopr*>   opr;          //Absvar as Functions of Stcvar
  
  int ivarOld;                       //Last abstract variable processed 
  
  OptActInfo** asTab;                //Array of abs-stc var. influence
  
 public:
  
  // watch VarType should correspond to type in STCVAR class
  // mpattr covers the case in which a variable affects attributes in
  // multiple physical domains
  
  
  enum VarTyp {notyp,attribute,coordinate,composite,nodalforce,
	       fluidparam,elecattr,thermattr,mixedvar,mpattr};
  
  
  Optvar();
  
  ~Optvar();
  
  void print(FILE*);
  void printres(FILE*);
  
  void buildAStable(int);

  void setvalptr(Structopt*);
  void setgradptr(Structopt*);
  
  void updvar();
  void resetvar();
  void initgrad();
  void updgrad(int ivar, double gradval);
  
  void setMeanValue(double* v=0);
  void generateRandom();
  void setFormFlag(int);
  
  double getFormL2norm2();
  double getgradFormL2norm2();
  
  double maxDesignVelocity(int iabs);

  size_t nAbs() const { return absvar.size(); }
  size_t nStc() const { return stcvar.size(); }
  int getIOldVar() const { return ivarOld; }
  Stcvar* getStcVar(size_t i) { return stcvar[i]; }
  Absvar* getAbsVar(size_t i) { return absvar[i]; }
  Optopr* getOpr(size_t i) { return opr[i]; }
  void addAbsVar(size_t, Absvar*);
  void addStcVar(size_t, Stcvar*);
  void addOpr(size_t, Optopr*);
};

//------------------------------------------------------------------------------

class OptActInfo {

     int lsize;

     int typ;

     ResizeArray<int> elements;

  public:
  		      
     OptActInfo();
  		      
     void add( int, int tp=Optvar::notyp ); 

     int check(int); 
     int getEntry(int);   

     int size()          { return lsize; }
     int getTyp()        { return typ; } 

     char* typname();

};

#endif

#endif
