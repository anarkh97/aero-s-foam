#ifndef _CONNECTIVITYT_H_
#define _CONNECTIVITYT_H_

#include <cstdio>
#include <map>
#include <iostream>
using std::map;

class BinFileHandler;

// component data structure
template<typename IndexType>
struct compStructT {
    IndexType numComp; // number of components
    IndexType *xcomp;  // pointer to renum for the beginning of each component
    IndexType *order;  // order of the nodes
    IndexType *renum;  // renumbering
    compStructT() { numComp=0; xcomp=0; order=0; renum=0; }
    void clearMemory() { if(xcomp) { delete [] xcomp; xcomp=0;}
        if(order) { delete [] order; order=0;}
        if(renum) { delete [] renum; renum=0;} }
};

//
// Access class used by BaseConnectivityT functions in order to free ConnectivityT from its link to element set
//

template <typename A>
struct base_traits;

template<typename A>
class DirectAccessT
{
    typedef typename base_traits<A>::IndexType IndexType;
    typedef typename base_traits<A>::DataType  DataType;

public :
    static IndexType getNum(A* oc, IndexType i)
    {return oc->num(i); }
    static DataType * getData(A* oc, IndexType i)
    {return (*oc)[i];}
    static IndexType getNumTarget(A *oc)
    {return oc->getNumTarget();}
    /*static IndexType * getTarget(A* oc)
	  {return oc->getTarget();}*/
    /*static IndexType * getPointer(A* oc)
	  {return oc->getPointer();}*/
    static IndexType getSize(A* oc)
    {return oc->csize();}
    // operator []
    /* static IndexType * hook(A* oc, IndexType i)
	   {return oc->operator[](i);}*/
};

template<typename IndexType, typename DataType> class ConnectivityT;

template<typename A, class Accessor = DirectAccessT<A> >
class BaseConnectivityT
{
    typedef typename base_traits<A>::IndexType IndexType;
    typedef typename base_traits<A>::DataType  DataType;

public :
    ConnectivityT<DataType,IndexType>* reverse(float *w = 0);
    template <class B, class AB>
    ConnectivityT<IndexType,typename base_traits<B>::DataType>* transcon(BaseConnectivityT<B,AB>* tc);
/*
  ConnectivityT<DataType,IndexType>* altReverse(float *w = 0);
  template <class B, class AB>
    ConnectivityT<IndexType,typename B::DataType>* altTranscon(BaseConnectivityT<B,AB>* tc);
  template <class B, class AB>
    ConnectivityT<IndexType,typename B::DataType>* transcon(BaseConnectivityT<B,AB>* tc, IndexType *, IndexType);
*/

    /* following functions are used to make transcon accept a BaseConnectivityT as argument */
    DataType *operator[](IndexType i) {return(Accessor::getData(static_cast<A*>(this),i));}
    IndexType csize() {return(Accessor::getSize(static_cast<A*>(this)));}
    IndexType num(IndexType i){return(Accessor::getNum(static_cast<A*>(this),i));}
    IndexType getNumTarget(){return(Accessor::getNumTarget(static_cast<A*>(this)));}
    void print()
    {
        for(IndexType i = 0; i < csize(); ++i )
        {
            std::cerr << i+1 << " -> ";
            IndexType numi = num(i);
            for(IndexType j = 0; j<numi; ++j)
                std::cerr << operator[](i)[j]+1 << " ";
            std::cerr << std::endl;
        }
    }
};

/** Class to access a set of element's individual connections */
template <typename A>
class SetAccessT
{
    typedef typename A::IndexType IndexType;
    typedef typename A::DataType  DataType;
public:
    IndexType size(); //<! returns the number of members of the set
    IndexType numNodes(IndexType i); //<! returns the number of targets for member i
    void getNodes(IndexType i, DataType *nd); //<! copies into nd the targets for member i
};

template<typename _IndexType = int, typename _DataType = int>
class ConnectivityT : public BaseConnectivityT<ConnectivityT<_IndexType,_DataType>,DirectAccessT<ConnectivityT<_IndexType,_DataType> > >
{
public:
    typedef _IndexType IndexType;
    typedef _DataType  DataType;

	IndexType getNumTarget() { return numtarget; }
    auto &getTarget()   { return target; }
    auto &getPointer() { return pointer; }

    ConnectivityT() { size = 0; numtarget = 0;}
    template <class A> ConnectivityT(SetAccessT<A> &sa);
    ConnectivityT(IndexType _size, IndexType *_pointer, DataType *_target, bool _removeable = true, float *_weight = 0);
    ConnectivityT(IndexType _size, IndexType *count);
    ConnectivityT(IndexType _size, IndexType count);
    ConnectivityT(BinFileHandler &, bool oldSower = false);
    ConnectivityT(IndexType ns); //dec
    ConnectivityT(FILE *f, IndexType nElem); // JAT 100614

    size_t write(BinFileHandler& f);
    size_t read(FILE* f);

    void countlink(IndexType from, IndexType to); //DEC
    void addlink(IndexType from, IndexType to); //DEC

    virtual ~ConnectivityT() = default;
    DataType *operator[](IndexType i);
    IndexType csize();
    IndexType numConnect(); // Total number of connections
    IndexType num(IndexType);
    IndexType offset(IndexType i) { return pointer[i]; } // Begining of ith part
    IndexType offset(IndexType i, DataType j); // returns a unique id for connection i to j
    IndexType cOffset(IndexType i, DataType j); // returns offset(i,j) - pointer[i]
    bool locate(IndexType i,IndexType j)
    { for(IndexType k=pointer[i]; k<pointer[i+1]; ++k) if(target[k] == j) return true; return false; }
    ConnectivityT<DataType,IndexType>* reverse(float *w = 0)
    { return(BaseConnectivityT<ConnectivityT>::reverse(w)); } // creates t->s from s->t
    DataType getTargetValue(IndexType i) { return target[i]; }
    auto &ptr() { return pointer; }
    auto &tgt() { return target; }
    bool isDiagonal(); // returns true if each target is only connected to itself

#if 0
    virtual void end_count();
        IndexType num(IndexType, IndexType*);
        ConnectivityT* transconOne(ConnectivityT*);
        void findPseudoDiam(IndexType *n1, IndexType *n2, IndexType *mask=0);
        IndexType  rootLS(IndexType root, IndexType *xls, IndexType *ls, IndexType &w, IndexType *mask=0);

        // Create a rooted level structure
        IndexType *renumSloan(IndexType *mask, IndexType &firstNum, IndexType *ren = 0);
        IndexType *renumRCM(IndexType *mask, IndexType &firstNum, IndexType *ren = 0);
        IndexType *renumSloan();
        compStructT renumByComponent(IndexType);
        void print(FILE * = stderr, IndexType node=-1);
        IndexType findMaxDist(IndexType *);
        ConnectivityT *merge(ConnectivityT *cn);
        // Create a copy of this connectivity without...
	ConnectivityT *collapse();
        ConnectivityT *subSection(bool *);
	ConnectivityT *trim(ConnectivityT *);
	ConnectivityT *copy();
	void sortTargets();
        void renumberTargets(IndexType *map);
        void renumberTargets(map<IndexType, IndexType> &);
        IndexType numNonZeroP();

        ConnectivityT *modify();
        ConnectivityT *modifyAlt();
        void combine(ConnectivityT *con2, IndexType *&cmap, IndexType *cmap2);  // adds con2 to this
        // adds all the entries in cmap (of size addSize)to each of the line in the current connectivity specified by entries in cmap
        // e.g. cmap = [2 3] and (*this)[2] = [1 2 4 5] (*this)[3] = [3 5], then (*this)[2] becomes [1 2 4 5 3]; (*this)[3] becomes [3 5 2]
        ConnectivityT *combineAll(IndexType addSize, IndexType *cmap); //ADDED FOR HEV PROBLEM, EC, 20070820 

        ConnectivityT * SCOTCH_graphPart(IndexType partnbr);*/
#endif
protected:
	IndexType size;           // size of pointer
	IndexType numtarget;      // size of target, number of Connections
	std::vector<IndexType> pointer;       // pointer to target
	std::vector<DataType> target;        // value of the connectivity
	std::vector<float> weight;      // weights of pointer (or 0)

};

template <typename _IndexType, typename _DataType>
struct base_traits<ConnectivityT<_IndexType,_DataType> > {
    typedef _IndexType IndexType;
    typedef _DataType DataType;
};

template<typename IndexType, typename DataType>
IndexType
ConnectivityT<IndexType,DataType>::csize() { return size; }

template<typename IndexType, typename DataType>
IndexType
ConnectivityT<IndexType,DataType>::num(IndexType n) { return (n < size) ? pointer[n+1] - pointer[n] : 0; }

template<typename IndexType, typename DataType>
DataType *
ConnectivityT<IndexType,DataType>::operator[](IndexType i) { return target.data() + pointer[i]; }

template<typename IndexType, typename DataType>
IndexType
ConnectivityT<IndexType,DataType>::numConnect() { return numtarget; }

template<typename IndexType, typename DataType>
bool
ConnectivityT<IndexType,DataType>::isDiagonal() { return (numtarget == size) ? true : false; }

/*
class CountedConnectivityT: public ConnectivityT {
     IndexType *cnt;
  public:
     CountedConnectivityT(IndexType ns);
     virtual ~CountedConnectivityT();
     IndexType count(IndexType n);
     void end_count();
     void remove(IndexType from, IndexType to);
};
*/

/** \breif  ImplicitConnectivity allows to use a type of data represented by class A
 * as a Connectivity without constructing the real ConnectivityT, but through clever use
 * of an appropriate Accessor.
 */
template <class A, class Accessor>
class ImplicitConnectivityT :
    public BaseConnectivityT<ImplicitConnectivityT<A,Accessor>/*, DirectAccessT<A>*/ >
{
public:
    typedef typename base_traits<Accessor>::IndexType IndexType;
    typedef typename base_traits<Accessor>::DataType  DataType;
private:
    A a;
    IndexType cacheVal, cacheSize;
    DataType *cacheTg;
    IndexType cacheNumVal, cacheNumIdx;
public:
    ImplicitConnectivityT(A _a) : a(_a) { cacheSize = 0; cacheTg = 0;
        cacheVal=-1; cacheNumIdx = -1; }
    ~ImplicitConnectivityT() { if(cacheTg) delete [] cacheTg; }
    IndexType csize() { return Accessor::getSize(*a); }
    IndexType num(IndexType i) {
        if(cacheNumIdx != i) {
            cacheNumVal = Accessor::getNum(*a,i);
            cacheNumIdx = i;
        }
        return cacheNumVal;
    }
    DataType *getNewCacheVal(IndexType j) {
        IndexType n = num(j);
        if(n > cacheSize) {
            IndexType newSize = std::max(n, 3*cacheSize/2);
            if(cacheTg) delete [] cacheTg;
            cacheTg = new DataType[newSize];
            cacheSize = newSize;
        }
        cacheVal = j;
        return(Accessor::getData(*a,j,cacheTg));
    }
    DataType *operator[](IndexType i) {
        if(cacheVal == i)
            return cacheTg;
        else
            return getNewCacheVal(i);
    }
    IndexType getNumTarget() {
        IndexType n = csize();
        IndexType res = 0;
        for(IndexType i = 0; i < n; ++i)
            res += num(i);
        return res;
    }
};

template <typename A, typename Accessor>
struct base_traits<ImplicitConnectivityT<A,Accessor> > {
    typedef typename base_traits<Accessor>::IndexType IndexType;
    typedef typename base_traits<Accessor>::DataType DataType;
};


#include <iostream>

// reverse() return a new connectivity that is the reverse of the present one
template<typename A, class Accessor>
ConnectivityT<typename base_traits<A>::DataType, typename base_traits<A>::IndexType> *
BaseConnectivityT<A,Accessor>::reverse(float * w)
{
    // The reverse connectivity has the same size as the original
    IndexType size = csize(); //Accessor::getSize(static_cast<A*>(this));
    IndexType numTarget = getNumTarget(); //Accessor::getNumTarget(static_cast<A*>(this));
    IndexType *res_target = new IndexType[numTarget];
    // Find the max of target
    DataType maxtarg = -1;
    for(IndexType i = 0; i < size; ++i) {
        IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
        DataType *tg = Accessor::getData(static_cast<A*>(this), i);
        for(IndexType j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
    }
    DataType res_size = maxtarg+1;
    DataType *res_pointer = new DataType[res_size+1];
    for(DataType i = 0; i <= res_size; ++i) res_pointer[i] = 0;

    // Now do a first pass to fill in res_pointer
    for(IndexType i=0; i < size; ++i) {
        IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
        DataType *tg = Accessor::getData(static_cast<A*>(this), i);
        for(IndexType j = 0; j < nTg; ++j) res_pointer[tg[j]]++;
    }

    for(DataType i = 1; i <= res_size; ++i) res_pointer[i] += res_pointer[i-1];

    // Second pass fills in target
    for(IndexType i = 0; i < size; ++i) {
        IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
        DataType *tg = Accessor::getData(static_cast<A*>(this), i);
        for(IndexType j = 0; j < nTg; ++j)
            res_target[--res_pointer[tg[j]]] = i;
    }

    ConnectivityT<typename A::DataType, typename A::IndexType> *res
        = new ConnectivityT<typename A::DataType, typename A::IndexType>(res_size, res_pointer, res_target, true, w);

    return res;
}

// Important NOTE: transcon cannot be called with the tc == this if tc is
// an Implicit ConnectivityT!!!
template<typename A, class Accessor>
template<class B, class AB>
ConnectivityT<typename base_traits<A>::IndexType,typename base_traits<B>::DataType>* BaseConnectivityT<A,Accessor>::transcon(BaseConnectivityT<B,AB>* tc)
{
    // First find the biggest target so we can size arrays correctly
    typename B::DataType tgmax=-1;

    IndexType size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
    for(typename B::IndexType i = 0; i < tc->csize(); ++i)
        for(typename B::IndexType j = 0; j < tc->num(i); ++j)
            tgmax = std::max(tgmax, (*tc)[i][j]);
    tgmax++; // Important adjustment

    // Now we can size the array that flags if a target has been visited
    IndexType *flags = new IndexType[tgmax];
    for(typename B::DataType i = 0; i < tgmax; ++i) flags[i] = -1;

    // Compute the new pointers
    IndexType *np = new IndexType[size+1];
    IndexType cp = 0;
    for(IndexType i = 0; i < size; ++i) {
        np[i] = cp;
        IndexType nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
        DataType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
        for(IndexType j = 0; j < nTg; ++j) {
            DataType intermed = tg[j];
            for(typename B::IndexType k = 0; k < tc->num(intermed); ++k)
                if(flags[(*tc)[intermed][k]] != i) {
                    flags[(*tc)[intermed][k]] = i;
                    cp++;
                }
        }
    }
    np[size] = cp;

    // Now allocate and fill the new target
    for(typename B::DataType i = 0; i < tgmax; ++i)
        flags[i] = -1;
    typename B::DataType *ntg = new typename B::DataType[cp];
    cp = 0;
    for(IndexType i = 0; i < size; ++i) {
        IndexType nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
        DataType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
        for(IndexType j = 0; j < nTg; ++j) {
            DataType intermed = tg[j];
            for (typename B::IndexType k = 0; k < tc->num(intermed); ++k)
                if(flags[(*tc)[intermed][k]] != i) {
                    flags[(*tc)[intermed][k]] = i;
                    ntg[cp] = (*tc)[intermed][k];
                    cp++;
                }
        }
    }
    delete [] flags;
    ConnectivityT<IndexType,typename B::DataType> *res = new ConnectivityT<IndexType,typename B::DataType>(size, np, ntg);
    return res;
}

#if 0
template<typename A, class Accessor>
template<class B, class AB>
ConnectivityT* BaseConnectivityT<A,Accessor>::transcon(BaseConnectivityT<B,AB>* tc, IndexType *ownedby, IndexType me)
{
 IndexType i,j,k;

 // First find the biggest target so we can size arrays correctly
 IndexType tgmax=-1;

 IndexType size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
 for(i = 0; i < tc->csize(); ++i) {
   if(ownedby[i] != me) continue;
   for(j = 0; j < tc->num(i); ++j)
     tgmax = std::max(tgmax, (*tc)[i][j]);
 }
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 IndexType *flags = new IndexType[tgmax];
 for(i = 0; i < tgmax; ++i) flags[i] = -1;

 // Compute the new pointers
 IndexType *np = new IndexType[size+1];
 IndexType cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   IndexType nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   IndexType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     IndexType intermed = tg[j];
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
 IndexType *ntg = new IndexType[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   IndexType nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   IndexType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     IndexType intermed = tg[j];
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
 ConnectivityT *res = new ConnectivityT(size, np, ntg);
 return res;
}

template<typename A, class Accessor>
ConnectivityT *
BaseConnectivityT<A,Accessor>::altReverse(float * w)
{
 // PJSA: 12-12-05 this version of reverse maintains the original ordering
 // The reverse connectivity has the same size as the original
 IndexType size = csize(); //Accessor::getSize(static_cast<A*>(this));
 IndexType numTarget = getNumTarget(); //Accessor::getNumTarget(static_cast<A*>(this));
 IndexType *res_target = new IndexType[numTarget];

 // Find the max of target
 IndexType maxtarg = 0;
 IndexType i;
 for(i=0; i < size; ++i) {
   IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
   IndexType *tg = Accessor::getData(static_cast<A*>(this), i);
   for(IndexType j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
 }
 IndexType res_size = maxtarg+1;
 IndexType *res_pointer = new IndexType[res_size+1];
 for(i = 0; i <= res_size; ++i) res_pointer[i] = 0;
                                                                                                                       
 // Now do a first pass to fill in res_pointer
 for(i=0; i < size; ++i) {
   IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
   IndexType *tg = Accessor::getData(static_cast<A*>(this), i);
   for(IndexType j = 0; j < nTg; ++j) res_pointer[tg[j]+1]++;
 }
                                                                                                                       
 IndexType *count = new IndexType[res_size];
 for(i = 1; i <= res_size; ++i) {
   res_pointer[i] += res_pointer[i-1];
   count[i-1] = 0;
 }

 // Second pass fills in target
 for(i=0; i < size; ++i) {
   IndexType nTg = Accessor::getNum(static_cast<A*>(this), i);
   IndexType *tg = Accessor::getData(static_cast<A*>(this), i);
   for(IndexType j = 0; j < nTg; ++j)
     res_target[res_pointer[tg[j]] + count[tg[j]]++] = i;
 }
 delete [] count;

 ConnectivityT *res = new ConnectivityT(res_size, res_pointer, res_target, 1, w);
 return res;
}

template<typename A, class Accessor>
template<class B, class AB>
ConnectivityT* BaseConnectivityT<A,Accessor>::altTranscon(BaseConnectivityT<B,AB>* tc)
{
 // PJSA 12-12-05 this version of transcon doesn't include connections with self unless entirely internal
 IndexType i,j,k;

 // First find the biggest target so we can size arrays correctly
 IndexType tgmax=-1;

 IndexType size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
 for(i = 0; i < tc->csize(); ++i)
   for(j = 0; j < tc->num(i); ++j)
     tgmax = std::max(tgmax, (*tc)[i][j]);
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 IndexType *flags = new IndexType[tgmax];
 for(i = 0; i < tgmax; ++i) flags[i] = -1;

 // Compute the new pointers
 IndexType *np = new IndexType[size+1];
 IndexType cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   IndexType nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   IndexType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     IndexType intermed = tg[j];
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
 IndexType *ntg = new IndexType[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   IndexType nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
   IndexType *tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
   for(j = 0; j < nTg; ++j) {
     IndexType intermed = tg[j];
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
 ConnectivityT *res = new ConnectivityT(size, np, ntg);
 return res;
}
#endif

template<typename IndexType, typename DataType>
template <class A>
ConnectivityT<IndexType,DataType>::ConnectivityT(SetAccessT<A> &sa)
{
    IndexType i;

    size = sa.size();

    // Find out the number of targets we will have
    pointer = new IndexType[size+1] ;
    IndexType pp = 0;
    for(i=0; i < size; ++i) {
        pointer[i] = pp;
        pp += sa.numNodes(i);
    }
    pointer[size] = pp;
    numtarget = pp;

    // Create the target array
    target = new DataType[pp];

    // Fill it in
    for(i=0; i < size; ++i) {
        sa.getNodes(i, target+pointer[i]);
    }
}

#ifdef _TEMPLATE_FIX_
#include <Utils.d/ConnectivityTImpl.h>
#endif

#endif
