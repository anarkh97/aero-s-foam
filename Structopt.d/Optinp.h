#ifndef _OPTINP_H_
#define _OPTINP_H_

#define MAXOPR 50

//------------------------------------------------------------------------------
struct nlpdata;
struct funcall;
struct funcdata;
struct graddata;

//------------------------------------------------------------------------------

struct funcall {
    int        typ;
    funcdata * fdata;
};
   
//------------------------------------------------------------------------------

struct funcdata {
    int        numopr;
    int        numgen;
    int*       numstat;
    int*       oprtyp;
    int*       oprnum;
    int*       gena;
    int*       genb;
    int*       gens;
    double*    a;
    double*    p;
    double*    b;
    funcall*   subfunc;
};

//------------------------------------------------------------------------------

struct dldata {

     int    node;
     int    typ;

     double refVal;
     double difVal;
};


//------------------------------------------------------------------------------

struct critdata {

    int     typ;

    int     i[5];
    double  d[5];

    int     ipsize;
    int*    ip;

    int      dlsize;
    dldata*  dl;

    void initialize() { typ=0; ipsize=0; ip=0; dlsize=0; dl=0; }
 
    void cleanup()    { if (ipsize) delete [] ip;
                        if (dlsize) delete [] dl; }
}; 

//------------------------------------------------------------------------------

struct criterion {

    int     num;
 
    critdata data;

    int     igen;
    int     gena;
    int     genb;
    int     gens;

    double  t;   
    int     anaId;

    void initialize() { num=0; igen=0; data.initialize(); t=-1; anaId=0; }

    void cleanup() { data.cleanup(); }

};    

//------------------------------------------------------------------------------

struct absvardata {
    int    num;
    int    igen;
    int    gena;
    int    genb;
    int    gens;
    int    dist;
    double val;
    double scl;
    double mean;
    double sdev;
    double low;
    double upp;
    double vin;
};


//------------------------------------------------------------------------------

struct solverdata {
    int        typ;
    nlpdata  * param;
    graddata * pgrad;
};

//------------------------------------------------------------------------------

struct nlpdata {
    int    ival [30];
    double rval [30];
    int    iflag[30];
    int    rflag[30];

    void initialize() 
    { int i; for (i=0;i<30;i++) {iflag[i]=0;rflag[i]=0; } }
};

//------------------------------------------------------------------------------

struct anadata {
    int    typ;
    int    ival[4];
    double rval[4];
};

//------------------------------------------------------------------------------

struct graddata {
    int    typ;
    int    mth;
    int    filter;
    int    epstyp;
    int    filterTyp;
    int    filterScale;
    int    numFilCrit;
    int    numGroups;
    int    numFilGrps;
    double epsval;
    double radius;
    double maxCount;
    double minExp;
    double maxExp;
    int*   filcritList; 
    int*   filGrpsList; 
};    


//------------------------------------------------------------------------------

struct funcdef {
    int        num;
    funcall    falldata;
}; 

//------------------------------------------------------------------------------

struct sumdata {
    int     oprtyp;
    int     oprnum;
    int     igen;
    int     gena;
    int     genb;
    int     gens;
    double  a;
    double  p;
    double  b;
};

//------------------------------------------------------------------------------

struct oprdata {
    int     oprtyp,oprnum;
};

//------------------------------------------------------------------------------

struct gendata {
    int a,b,s;
}; 

//------------------------------------------------------------------------------

struct stcelvdata {
    int var,ele1,ele2,typ;
}; 

//------------------------------------------------------------------------------

struct intlist {

    int  mxsize;
    int  ipsize;
    int* ip;

    void add(int ival) {
      if (ipsize == mxsize) {
        mxsize = mxsize +10;
        int * newip=new int[mxsize];
        int j;
        for (j=0;j<ipsize;j++) newip[j]=ip[j];
        delete [] ip;
        ip = newip;
      }
      ip[ipsize]=ival;
      ipsize++;
    }    
        
    void cleanup() { delete [] ip; }
}; 

//------------------------------------------------------------------------------

struct dllist {

    int     mxsize;
    int     dlsize;
    dldata* dl;

    void add(int n, int t, double & dv, double & rv) 
    {
      if (dlsize == mxsize) {
        mxsize = mxsize +10;
        dldata * newdl=new dldata[mxsize];
        int j;
        for (j=0;j<dlsize;j++) {
           newdl[j]=dl[j]; }
        delete [] dl;
        dl = newdl;
      }
      dl[dlsize].node   = n;
      dl[dlsize].typ    = t;
      dl[dlsize].refVal = rv;
      dl[dlsize].difVal = dv;

      dlsize++;
    }    
        
    void cleanup() { delete [] dl; }

};

//------------------------------------------------------------------------------

extern funcdata* buildFuncdata();

extern void initFuncall  (funcall&);
extern void cleanFuncall (funcall&);
extern void cleanFuncdata(funcdata&);

#endif
