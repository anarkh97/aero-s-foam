#ifndef _STRUCTPROP_H_
#define _STRUCTPROP_H_

#include <Element.d/Element.h>
#include <Driver.d/Mpc.h>
#include <map>

//-------------------------------------------------------------

// Element attribute structure
struct Attrib {
  int nele;
  int attr;
  int cmp_attr, cmp_frm;
  double cmp_theta;  // PJSA 4-8-05: angle between element edge01 and material x axis
  Attrib() { nele = 0; attr = cmp_attr = cmp_frm = -1; cmp_theta = 0.0; }
};

class ElemAttrAccessor{
 public:
  static int getNum(/*const*/ std::pair<int,std::map<int,Attrib> *> &, int )
    { return 1; }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
     { 
       if(nd) { nd[0] = (*(o.second))[i].nele; return nd; }
       else return &(*(o.second))[i].nele;
     }
};

class MatAttrAccessor{ // EleToMat
 public:
  static int getNum(/*const*/ std::pair<int,std::map<int,Attrib> *> &, int )
    { return 1; }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
     { 
       if(nd) { nd[0] = (*(o.second))[i].attr; return nd; }
       else return &(*(o.second))[i].attr;
     }
};



class CmpAttrAccessor{
 public:
  static int getNum(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i)
    { 
      if ( (*(o.second))[i].cmp_attr != -1 )
	return 1; 
      else
	return 0;
    }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
     { 
       if ( (*(o.second))[i].cmp_attr != -1 )
	 {
	   if(nd) { nd[0] = (*(o.second))[i].cmp_attr; return nd; }
	   else return &(*(o.second))[i].cmp_attr;
	 }
       std::cerr << "warning getting data for size = 0 !" << std::endl;
       return 0;
     }
};

class CmpFrAttrAccessor{
 public:
  static int getNum(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i)
    { 
      if ( (*(o.second))[i].cmp_frm != -1 )
	return 1; 
      else
	return 0;
    }
  static int getSize(const std::pair<int,std::map<int,Attrib> *> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,std::map<int,Attrib> *> &o, int i, int *nd)
     { 
       if ( (*(o.second))[i].cmp_frm != -1 )
	 {
	   if(nd) { nd[0] = (*(o.second))[i].cmp_frm; return nd; }
	   else return &(*(o.second))[i].cmp_frm;
	 }
       std::cerr << "warning getting data for size = 0 !" << std::endl;
       return 0;
     }
};

class LMPCAccessor{
 public:
  static int getNum(/*const*/ std::pair<int,ResizeArray<LMPCons*> *> &o, int i)
    { 
      if ( (*(o.second))[i]->nterms > 0 )
	return (*(o.second))[i]->nterms; 
      else
	return 0;
    }
  static int getSize(const std::pair<int,ResizeArray<LMPCons*> *> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,ResizeArray<LMPCons*> *> &o, int i, int *nd)
     { 
       if ( (*(o.second))[i]->nterms > 0 )
	 {
	   if(nd) 
	     {
	       for(int j = 0 ; j<(*(o.second))[i]->nterms ; ++j) // potential memory leak
		 nd[j] = (*(o.second))[i]->terms[j].nnum;
	       return nd; 
	     }
	   else
	     {
	       nd = new int[(*(o.second))[i]->nterms];
	       for(int j = 0 ; j<(*(o.second))[i]->nterms ; ++j)
		 nd[j] = (*(o.second))[i]->terms[j].nnum;
	       return nd;
	     }
	 }
       std::cerr << "warning getting data for size = 0 !" << std::endl;
       return 0;
     }
};


#include <Utils.d/BlockAlloc.h>
typedef std::map<int, StructProp, std::less<int>, block_allocator<std::pair<const int, StructProp> > > SPropContainer;


// NOT THREAD SAFE...
// BUT SOWER IS ON ONE THREAD
class CurveMatAccessor{
  // protect StructProp->__SOWER_TMP for multi-threading TG
 public:
  static int getNum(/*const*/ std::pair<int, SPropContainer*> &o, int i)
    { 
      if ( ((int) - (*(o.second))[i].W) > 0 )
	return 1; 
      else
	return 0;
    }
  static int getSize(const std::pair<int, SPropContainer*> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int, SPropContainer*> &o, int i, int *nd)
     { 
       if ( ((int) - (*(o.second))[i].W) > 0 )
	 {
	   if(nd) 
	     { 
	       nd[0] = (int) - (*(o.second))[i].W; 
	       return nd; 
	     }
	   else 
	     {
	       (*(o.second))[i].__SOWER_TMP = (int) - (*(o.second))[i].W;
	       return &((*(o.second))[i].__SOWER_TMP);
	     }
	 }
       std::cerr << "warning getting data for size = 0 !" << std::endl;
       return 0;
     }
};

// NOT THREAD SAFE...
// BUT SOWER IS ON ONE THREAD
class CurveYoungMatAccessor{
  // protect StructProp->__SOWER_TMP for multi-threading TG
 public:
  static int getNum(/*const*/ std::pair<int, SPropContainer *> &o, int i)
    { 
      if ( ((int) - (*(o.second))[i].E) > 0 )
	return 1; 
      else
	return 0;
    }
  static int getSize(const std::pair<int,SPropContainer*> &o)
    { return o.first; }
  static int *getData(/*const*/ std::pair<int,SPropContainer*> &o, int i, int *nd)
     { 
       if ( ((int) - (*(o.second))[i].E) > 0 )
	 {
	   if(nd) 
	     { 
	       nd[0] = (int) - (*(o.second))[i].E; 
	       return nd; 
	     }
	   else 
	     {
	       (*(o.second))[i].__SOWER_TMP = (int) - (*(o.second))[i].E;
	       return &((*(o.second))[i].__SOWER_TMP);
	     }
	 }
       std::cerr << "warning getting data for size = 0 !" << std::endl;
       return 0;
     }
};



//-------------------------------------------------------------

#include <Driver.d/EFrameData.h> // TG : struct EFrameData moved here

//-------------------------------------------------------------

struct ControlInfo {
   FILE *checkfileptr;
   char *checkfile;
   char *nodeSetName;
   char *elemSetName;
   char *bcondSetName;
   char *decomposition;
   FILE *decPtr;
   char *currentRestartFile;
   char *lastRestartFile;
   char *outputExt;
   char *FlagRST;
   ControlInfo() { checkfile   = (char *) "check"; nodeSetName   = (char *) "nodeset";
                   elemSetName = (char *) "elemset"; bcondSetName = (char *) "bcondset";
                   decomposition = (char *) "DECOMPOSITION";
                   decPtr = 0; checkfileptr = 0;
                   currentRestartFile = 0;
                   lastRestartFile    = 0; outputExt = (char *) ""; FlagRST = (char *) "old"; }
};

#endif
