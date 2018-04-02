#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

#include <cstdio>
#include <map>
#include <iostream>

class Elemset;
class EqNumberer;
class BinFileHandler;
class SommerElement;
class FaceElemSet;

// component data structure
struct compStruct {
  int numComp; // number of components
  int *xcomp;  // pointer to renum for the beginning of each component 
  int *order;  // order of the nodes
  int *renum;  // renumbering
  compStruct() { numComp = 0; xcomp=0; order=0; renum=0; }
  void clearMemory() { if(xcomp) { delete [] xcomp; xcomp=0;}
                       if(order) { delete [] order; order=0;}
                       if(renum) { delete [] renum; renum=0;} }
};

/*
 * Access class used by BaseConnectivity functions in order to free Connectivity from its link to element set
*/
template<typename A>
class DirectAccess {
  public : 
  static int getNum(A* oc, int i)
    {return oc->num(i); }
  static int * getData(A* oc, int i)
    {return (*oc)[i];}
  static int getNumTarget(A *oc)
    {return oc->getNumTarget();}
  /*static int * getTarget(A* oc)
    {return oc->getTarget();}*/
  /*static int * getPointer(A* oc)
    {return oc->getPointer();}*/
 static int getSize(A* oc)
   {return oc->csize();}
 // operator []
 /* static int * hook(A* oc, int i)
    {return oc->operator[](i);}*/
};

class Connectivity;

template<typename A, class Accessor = DirectAccess<A> >
class BaseConnectivity
{
  public :
  Connectivity* reverse(float *w = 0);
  Connectivity* altReverse(float *w = 0);
  template <class B, class AB>
    Connectivity* transcon(BaseConnectivity<B,AB>* tc);
  template <class B, class AB>
    Connectivity* altTranscon(BaseConnectivity<B,AB>* tc);
  template <class B, class AB>
    Connectivity* transcon(BaseConnectivity<B,AB>* tc, int *, int);
  
  /* following functions are used to make transcon accept a BaseConnectivity as argument */
  int *operator[](int i) {return(Accessor::getData(static_cast<A*>(this),i));}
  int csize() {return(Accessor::getSize(static_cast<A*>(this)));}
  int num(int i){return(Accessor::getNum(static_cast<A*>(this),i));}
  int getNumTarget(){return(Accessor::getNumTarget(static_cast<A*>(this)));}
  void print()
  {
    for(int i = 0; i < csize(); ++i )
    {
      std::cerr << i+1 << " -> " ; 
      int numi = num(i);
      for(int j = 0; j<numi ; ++j)
      std::cerr << operator[](i)[j] << ",";
      std::cerr << std::endl;
    }
  }
};

/** Class to access a set of element's individual connections */
template <class A>
class SetAccess {
public:
	int size(); //<! returns the number of members of the set
	int numNodes(int i); //<! returns the number of targets for member i
	void getNodes(int i, int *nd); //<! copies into nd the targets for member i
};

class Connectivity : public BaseConnectivity<Connectivity,DirectAccess<Connectivity> >
{
 protected:
        int removeable;     // whether to delete pointer memory or not
	int size;           // size of pointer
	int numtarget;      // size of target, number of Connections
	int *pointer;       // pointer to target
	int *target;        // value of the connectivity
	float *weight;      // weights of pointer (or 0)
       
 public:
        typedef int IndexType;
	int getNumTarget() {return numtarget; }
	int * getTarget() {return target; }
	int * getPointer() {return pointer; }
	/*	void displayRelationships(void)
	  {
	    for(int i = 0; i< size; i++)
	      {
		int bot = pointer[i];
		int top = pointer[i+1]-bot;
		cout << "element " << i << " voisins : ";
		for(int j = bot; j < top; j++)
		  cout << target[j] << "," ;
		cout << endl;
	      }
 	    
	      }*/
        Connectivity() { size = 0; numtarget = 0; pointer = 0; target = 0; weight = 0; removeable = 1; }
        template <class A> Connectivity(SetAccess<A> &sa);
        Connectivity(Elemset *);
        Connectivity(Elemset *, int, SommerElement**);
        Connectivity(int _size, int *_pointer, int *_target, int _removeable=1, float *_weight = 0);
        Connectivity(int _size, int *count);
        Connectivity(int _size, int count);
	Connectivity(BinFileHandler &, bool oldSower = false);
        Connectivity(FaceElemSet*, int size = 0);
	Connectivity(int ns); //dec
	Connectivity(FILE *f, int nElem); // JAT 100614
        Connectivity(Elemset *els, Connectivity *nodeToElem);
        Connectivity(const Connectivity&);
	size_t write(BinFileHandler& f);
	size_t writeg(BinFileHandler& f);
	size_t read(FILE* f);
        size_t write(FILE* f);

  	void countlink(int from, int to); //DEC
        void addlink(int from, int to); //DEC

        Connectivity(SommerElement  **, int);
        Connectivity(SommerElement  **, int, int);

        virtual ~Connectivity();
	virtual void end_count(); //dec
        void setRemoveable(int _removeable) { removeable = _removeable; }
        int *operator[](int i);
        int csize();
        int numConnect(); // Total number of connections
        int num(int);
        int num(int, int*);
        int offset(int i) { return pointer[i]; } // Begining of ith part
        int offset(int i,int j); // returns a unique id for connection i to j
        int cOffset(int i,int j); // returns offset(i,j) - pointer[i]
        bool locate(int i,int j) 
          { for(int k=pointer[i]; k<pointer[i+1]; ++k) if(target[k] == j) return true; return false; }
	// this call to ::reverse makes HB's mortar lib compile... investigate please.
        Connectivity* reverse(float *w = 0)
          { return(BaseConnectivity<Connectivity>::reverse(w)); } // creates t->s from s->t
        Connectivity* transconOne(Connectivity*);
        int getTargetValue(int i) { return target[i]; }

        void findPseudoDiam(int *n1, int *n2, int *mask=0);
        int  rootLS(int root, int *xls, int *ls, int &w, int *mask=0);

long long memsize() {return ((long long)size + numtarget + 1)*sizeof(int);} 

        // Create a rooted level structure
        int *renumSloan(int *mask, int &firstNum, int *ren = 0);
        int *renumRCM(int *mask, int &firstNum, int *ren = 0);
        int *renumSloan();
        compStruct renumByComponent(int);
        void print(FILE * = stderr, int node=-1);
        int findMaxDist(int *);
        int findProfileSize(EqNumberer *eqNumber, int unroll=1);
        int findProfileSizes(EqNumberer *eqNumber, compStruct &, long *,
                             int unroll=1);
        int *ptr() { return pointer; }
        int *tgt() { return target; }
        Connectivity *merge(Connectivity *cn);
        // Create a copy of this connectivity without...
	Connectivity *collapse();
        Connectivity *subSection(bool *);
	Connectivity *trim(Connectivity *);
	Connectivity *copy();
	void sortTargets();
        void renumberTargets(int *map);
        void renumberTargets(std::map<int, int> &);
        int numNonZeroP();

        bool isDiagonal(); // returns true if each target is only connected to itself
        Connectivity *modify();
        Connectivity *modifyAlt();
        void combine(Connectivity *con2, int *&cmap, int *cmap2);  // adds con2 to this
        // adds all the entries in cmap (of size addSize)to each of the line in the current connectivity specified by entries in cmap
        // e.g. cmap = [2 3] and (*this)[2] = [1 2 4 5] (*this)[3] = [3 5], then (*this)[2] becomes [1 2 4 5 3]; (*this)[3] becomes [3 5 2]
        Connectivity *combineAll(int addSize, int *cmap);

        double estimateComponentCost(EqNumberer *eqn, compStruct &cs, double *cost, double *bandwidth, 
                                     double coef=400, int unroll=1);
        double estimateCost(EqNumberer *eqn, double &cost, double &bandwidth,
                            double coef=400, int unroll=1);

        Connectivity * SCOTCH_graphPart(int partnbr);
};



inline int
Connectivity::csize() { return size; }

inline int
Connectivity::num(int n) { return (n < size) ? pointer[n+1] - pointer[n] : 0; }

inline int *
Connectivity::operator[](int i) { return target +pointer[i] ; }

inline int
Connectivity::numConnect() { return numtarget; }

inline bool
Connectivity::isDiagonal() { return (numtarget == size) ? true : false; }

class CountedConnectivity: public Connectivity {
     int *cnt;
  public:
     CountedConnectivity(int ns);
     virtual ~CountedConnectivity();
     int count(int n);
     void end_count();
     void remove(int from, int to);
};

/*
 * class Implicit connectivity will allow to use a type of data represented by class A
 * as a Connectivity without constructing the real Connectivity, but through clever use
 * of an appropriate Accessor
 */
template <class A, class Accessor>
class ImplicitConnectivity :
public BaseConnectivity<ImplicitConnectivity<A,Accessor>/*, DirectAccess<A>*/ > {
    A a;
    int cacheVal, cacheSize;
    int *cacheTg;
    int cacheNumVal, cacheNumIdx;
  public:
    ImplicitConnectivity(A _a) : a(_a) { cacheSize = 0; cacheTg = 0;
                                         cacheVal=-1; cacheNumIdx = -1; }
    ~ImplicitConnectivity() { if(cacheTg) delete [] cacheTg; }
    int csize() { return Accessor::getSize(*a); }
    int num(int i) { if(cacheNumIdx != i) {
                        cacheNumVal = Accessor::getNum(*a,i);
                        cacheNumIdx = i;
                     }
                     return cacheNumVal;
                   }
    int *getNewCacheVal(int j) {
      int n = num(j);
      if(n > cacheSize) {
         int newSize = std::max(n, 3*cacheSize/2);
	 if(cacheTg) delete [] cacheTg;
	 cacheTg = new int[newSize];
	 cacheSize = newSize;
      }
      cacheVal = j;
      return(Accessor::getData(*a,j,cacheTg));
    }
    int *operator[](int i) {
      if(cacheVal == i)
        return cacheTg;
      else
        return getNewCacheVal(i);
    }
    int getNumTarget() {
       int n = csize();
       int res = 0;
       for(int i = 0; i < n; ++i)
         res += num(i);
       return res;
    }
};

#include <iostream>

// reverse() return a new connectivity that is the reverse of the present one
template<typename A, class Accessor>
Connectivity * 
BaseConnectivity<A,Accessor>::reverse(float * w)
{
 // The reverse connectivity has the same size as the original
 int size = csize(); //Accessor::getSize(static_cast<A*>(this));
 int numTarget = getNumTarget(); //Accessor::getNumTarget(static_cast<A*>(this));
 int *res_target = new int[numTarget];

 // Find the max of target
 int maxtarg = -1; // PJSA
 int i;
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
 }
 int res_size = maxtarg+1;
 int *res_pointer = new int[res_size+1];
 for(i = 0; i <= res_size; ++i) res_pointer[i] = 0;

 // Now do a first pass to fill in res_pointer
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j) res_pointer[tg[j]]++;
 }

 for(i = 1; i <= res_size; ++i) res_pointer[i] += res_pointer[i-1];

 // Second pass fills in target
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j)
     res_target[--res_pointer[tg[j]]] = i;
 }

 Connectivity *res = new Connectivity(res_size, res_pointer, res_target, 1, w);
 return res;
}

// Important NOTE: transcon cannot be called with the tc == this if tc is
// an Implicit Connectivity!!!
template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::transcon(BaseConnectivity<B,AB>* tc)
{
 int i,j,k;

 // First find the biggest target so we can size arrays correctly
 int tgmax=-1;

 int size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
 for(i = 0; i < tc->csize(); ++i)
   for(j = 0; j < tc->num(i); ++j)
     tgmax = std::max(tgmax, (*tc)[i][j]);
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 int *flags = new int[tgmax];
 for(i = 0; i < tgmax; ++i) flags[i] = -1;
 
 // Compute the new pointers
 int *np = new int[size+1];
 int cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     for(k = 0; k < tc->num(intermed); ++k)
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         cp ++;
       }
   }
 }
 np[size] = cp;

 // Now allocate and fill the new target
 for(i = 0; i < tgmax; ++i)
   flags[i] = -1;
 int *ntg = new int[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     for (k = 0; k < tc->num(intermed); ++k)
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         ntg[cp] = (*tc)[intermed][k];
         cp ++;
       }
   }
 }
 delete [] flags;
 Connectivity *res = new Connectivity(size, np, ntg);
 return res;
}

template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::transcon(BaseConnectivity<B,AB>* tc, int *ownedby, int me)
{
 int i,j,k;

 // First find the biggest target so we can size arrays correctly
 int tgmax=-1;

 int size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
 for(i = 0; i < tc->csize(); ++i) {
   if(ownedby[i] != me) continue;
   for(j = 0; j < tc->num(i); ++j)
     tgmax = std::max(tgmax, (*tc)[i][j]);
 }
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 int *flags = new int[tgmax];
 for(i = 0; i < tgmax; ++i) flags[i] = -1;

 // Compute the new pointers
 int *np = new int[size+1];
 int cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     if(ownedby[intermed] != me) continue;
     for(k = 0; k < tc->num(intermed); ++k)
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         cp ++;
       }
   }
 }
 np[size] = cp;

 // Now allocate and fill the new target
 for(i = 0; i < tgmax; ++i)
   flags[i] = -1;
 int *ntg = new int[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     if(ownedby[intermed] != me) continue;
     for (k = 0; k < tc->num(intermed); ++k)
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         ntg[cp] = (*tc)[intermed][k];
         cp ++;
       }
   }
 }
 delete [] flags;
 Connectivity *res = new Connectivity(size, np, ntg);
 return res;
}

template<typename A, class Accessor>
Connectivity *
BaseConnectivity<A,Accessor>::altReverse(float * w)
{
 // PJSA: 12-12-05 this version of reverse maintains the original ordering
 // The reverse connectivity has the same size as the original
 int size = csize(); //Accessor::getSize(static_cast<A*>(this));
 int numTarget = getNumTarget(); //Accessor::getNumTarget(static_cast<A*>(this));
 int *res_target = new int[numTarget];

 // Find the max of target
 int maxtarg = 0;
 int i;
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
 }
 int res_size = maxtarg+1;
 int *res_pointer = new int[res_size+1];
 for(i = 0; i <= res_size; ++i) res_pointer[i] = 0;
                                                                                                                       
 // Now do a first pass to fill in res_pointer
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j) res_pointer[tg[j]+1]++;
 }
                                                                                                                       
 int *count = new int[res_size];
 for(i = 1; i <= res_size; ++i) {
   res_pointer[i] += res_pointer[i-1];
   count[i-1] = 0;
 }

 // Second pass fills in target
 for(i=0; i < size; ++i) {
   int nTg = Accessor::getNum(static_cast<A*>(this), i);
   int *tg = Accessor::getData(static_cast<A*>(this), i);
   for(int j = 0; j < nTg; ++j)
     res_target[res_pointer[tg[j]] + count[tg[j]]++] = i;
 }
 delete [] count;

 Connectivity *res = new Connectivity(res_size, res_pointer, res_target, 1, w);
 return res;
}

template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::altTranscon(BaseConnectivity<B,AB>* tc)
{
 // PJSA 12-12-05 this version of transcon doesn't include connections with self unless entirely internal
 int i,j,k;
                                                                                                                       
 // First find the biggest target so we can size arrays correctly
 int tgmax=-1;
                                                                                                                       
 int size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
 for(i = 0; i < tc->csize(); ++i)
   for(j = 0; j < tc->num(i); ++j)
     tgmax = std::max(tgmax, (*tc)[i][j]);
 tgmax++; // Important adjustment
                                                                                                                       
 // Now we can size the array that flags if a target has been visited
 int *flags = new int[tgmax];
 for(i = 0; i < tgmax; ++i) flags[i] = -1;
                                                                                                                       
 // Compute the new pointers
 int *np = new int[size+1];
 int cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     for(k = 0; k < tc->num(intermed); ++k) {
       if(((*tc)[intermed][k] == i) && (tc->num(intermed) != 1)) continue;
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         cp ++;
       }
     }
   }
 }
 np[size] = cp;
                                                                                                                       
 // Now allocate and fill the new target
 for(i = 0; i < tgmax; ++i)
   flags[i] = -1;
 int *ntg = new int[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   int *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     int intermed = tg[j];
     for (k = 0; k < tc->num(intermed); ++k) {
       if(((*tc)[intermed][k] == i) && (tc->num(intermed) != 1)) continue;
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         ntg[cp] = (*tc)[intermed][k];
         cp ++;
       }
     }
   }
 }
 delete [] flags;
 Connectivity *res = new Connectivity(size, np, ntg);
 return res;
}

template <class A>
Connectivity::Connectivity(SetAccess<A> &sa)
{
 removeable = 1;
 int i;
 weight = (float *) 0;

 size = sa.size();

 // Find out the number of targets we will have
 pointer = new int[size+1] ;
 int pp = 0;
 for(i=0; i < size; ++i) {
   pointer[i] = pp;
   pp += sa.numNodes(i);
 }
 pointer[size] = pp;
 numtarget = pp;

 // Create the target array
 target = new int[pp];

 // Fill it in
 for(i=0; i < size; ++i) {
   sa.getNodes(i, target+pointer[i]);
 }
}
#endif
