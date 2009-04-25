#ifndef _OPTOPR_H_
#define _OPTOPR_H_

#ifdef STRUCTOPT

#include <stdio.h>

#include <Utils.d/resize_array.h>

class Absvar;
class Optcrit;
class OptActInfo;

struct funcall;

//------------------------------------------------------------------------------
//                      Functions on Optimization Criteria
//                                                       created  9/1/98 by Kurt
//------------------------------------------------------------------------------
class Optopr 
{  
 public:  
  //  Operator :          0 - Weighted Sum of Criteria
  //                      1 - Kreiselmeier-Steinhauser
  //                      2 - Weighted Multiplication
  //                      3 - SINUS
  //                      4 - COSINUS
  //                      5 - Exponential function
  //                      6 - User defined function  
  int typ;
 public:
  virtual ~Optopr() {}
  char * gettypname();
  
  static Optopr * create( int _typ=0 );
  
  virtual void build ( Optopr  *_opr, double a, double b, double c)=0;
  virtual void build ( Optcrit *_opc, double a, double b, double c)=0;
  virtual void build ( Absvar  *_opv, double a, double b, double c)=0;
  
  virtual void print    (FILE*)=0;	
  virtual void printres (FILE*)=0;	
  
  virtual double getval  ()=0;
  virtual double getgrad ()=0;
  
  virtual void extractAbsvar(OptActInfo&)=0;
  virtual void extractCriteria(int *)=0;
  
  virtual funcall & extractFuncall()=0; 
};

//------------------------------------------------------------------------------
//      Operator: sum ( a[i] * criteria [i]^ p[i] + b [i] )
//------------------------------------------------------------------------------

class OptoprSum : public Optopr 
{ 
 public:  
  int numsum;                             //Number of Summands  
  int numopc;                             //Number of Criteria
  int numopr;                             //Number of Functions
  int numopv;                             //Number of Variables
  
  ResizeArray<Optcrit*>  opc;             //Criteria
  ResizeArray<Optopr* >   opr;             //Functions
  ResizeArray<Absvar*>   opv;             //Variables
  
  ResizeArray<int>   indopc;              //Criteria
  ResizeArray<int>   indopr;              //Functions
  ResizeArray<int>   indopv;              //Variables
  
  ResizeArray<double>    a;               //Factor
  ResizeArray<double>    b;               //Constant
  ResizeArray<double>    p;               //Exponent
  
  
 public:
  
  OptoprSum ();
  
  void print    (FILE*);
  void printres (FILE*);
  
  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
  
  double getval  ();
  double getgrad (); 
  
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);
  
  funcall & extractFuncall();
};

//------------------------------------------------------------------------------
//      Operator: ln sum ( e ^ ( a[i] * criteria [i] ^ p[i] + b [i] ) )
//------------------------------------------------------------------------------
class OptoprKS : public Optopr 
{
 public: 
  int numsum;                             //Number of Summands
  
  int numopc;                             //Number of Criteria
  int numopr;                             //Number of Functions
  int numopv;                             //Number of Variables
  
  ResizeArray<Optcrit*>  opc;             //Criteria
  ResizeArray<Optopr* >   opr;             //Functions
  ResizeArray<Absvar*>   opv;             //Variables
  
  ResizeArray<int>   indopc;              //Criteria
  ResizeArray<int>   indopr;              //Functions
  ResizeArray<int>   indopv;              //Variables
  
  ResizeArray<double>    a;               //Factor
  ResizeArray<double>    b;               //Constant
  ResizeArray<double>    p;               //Exponent
  
  
 public:
  
  OptoprKS ();
  
  void print    (FILE*);
  void printres (FILE*);
  
  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
  
  double getval  ();
  double getgrad (); 
  
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);
  
  funcall & extractFuncall();
};

//------------------------------------------------------------------------------
//      Operator: PI ( a[i] * criteria [i] ^ p[i] + b [i] )
//------------------------------------------------------------------------------
class OptoprMul : public Optopr 
{  
 public:
  
  int numsum;                             //Number of Summands
  
  int numopc;                             //Number of Criteria
  int numopr;                             //Number of Functions
  int numopv;                             //Number of Variables
  
  ResizeArray<Optcrit*>  opc;             //Criteria
  ResizeArray<Optopr* >   opr;             //Functions
  ResizeArray<Absvar*>   opv;             //Variables
  
  ResizeArray<int>   indopc;              //Criteria
  ResizeArray<int>   indopr;              //Functions
  ResizeArray<int>   indopv;              //Variables
  
  ResizeArray<double>    a;               //Factor
  ResizeArray<double>    b;               //Constant
  ResizeArray<double>    p;               //Exponent
  
  
 public:
  
  OptoprMul ();
  
  void print    (FILE*);
  void printres (FILE*);

  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
     
  double getval  ();
  double getgrad (); 
       
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);

  funcall & extractFuncall();
};

//------------------------------------------------------------------------------
//      Operator:  sin( a * criteria ^ p + b ) ; cos( a * criteria ^ p + b )
//------------------------------------------------------------------------------
class OptoprTrig : public Optopr 
{ 
 public:
          
  Optcrit*  opc;             //Criteria
  Optopr*   opr;             //Functions
  Absvar*   opv;             //Variables

  double    a;               //Factor
  double    b;               //Constant
  double    p;               //Exponent
       
      
 public:

  OptoprTrig ();

  void print    (FILE*);
  void printres (FILE*);

  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
     
  double getval  ();
  double getgrad (); 
       
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);

  funcall & extractFuncall();
};

//------------------------------------------------------------------------------
//      Operator:   MUL  of exp( a * criteria ^ p + b ) 
//------------------------------------------------------------------------------
class OptoprExp : public Optopr 
{
 public:
 
  int numsum;                             //Number of Summands
 
  int numopc;                             //Number of Criteria
  int numopr;                             //Number of Functions
  int numopv;                             //Number of Variables
          
  ResizeArray<Optcrit*>  opc;             //Criteria
  ResizeArray<Optopr* >   opr;             //Functions
  ResizeArray<Absvar*>   opv;             //Variables

  ResizeArray<int>   indopc;              //Criteria
  ResizeArray<int>   indopr;              //Functions
  ResizeArray<int>   indopv;              //Variables

  ResizeArray<double>    a;               //Factor
  ResizeArray<double>    b;               //Constant
  ResizeArray<double>    p;               //Exponent
       
      
 public:

  OptoprExp ();

  void print    (FILE*);
  void printres (FILE*);

  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
     
  double getval  ();
  double getgrad (); 
       
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);

  funcall & extractFuncall();
};

//------------------------------------------------------------------------------
//      Operator:   User-defined function 
//------------------------------------------------------------------------------
class OptoprUdef : public Optopr 
{
 public:
     
  int idUdef;                             //Id-number of user function

  int numsum;                             //Number of Summands
 
  int numopc;                             //Number of Criteria
  int numopr;                             //Number of Functions
  int numopv;                             //Number of Variables
          
  ResizeArray<Optcrit*>  opc;             //Criteria
  ResizeArray<Optopr* >   opr;             //Functions
  ResizeArray<Absvar*>   opv;             //Variables

  ResizeArray<int>   indopc;              //Criteria
  ResizeArray<int>   indopr;              //Functions
  ResizeArray<int>   indopv;              //Variables

  ResizeArray<int>   indpar;              //Position of in paramter list
       
  void * handle;                          //Handle to share object

  double (*fct)(double*,double*,int,int); //Pointer to user function
      
 public:

  OptoprUdef ();

  void print    (FILE*);
  void printres (FILE*);

  void build ( Optopr  *_opr, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Optcrit *_opc, double _a=1.0, double _p=1.0, double _b=0.0);
  void build ( Absvar  *_opv, double _a=1.0, double _p=1.0, double _b=0.0);
     
  double getval  ();
  double getgrad (); 
       
  void extractAbsvar(OptActInfo&);
  void extractCriteria(int *);

  funcall & extractFuncall();
       
  void getUdef();
};

#endif
 
#endif
