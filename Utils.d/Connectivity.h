#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

#include <cstdio>
#include <map>
#include <vector>
#include <iostream>
#include <Extermal.d/include/gsl/span>

class Elemset;
class EqNumberer;
class BinFileHandler;
class SommerElement;

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

inline
const int *__getPtr(const int *p) { return p; }

template <typename T>
const T *__getPtr(const gsl::span<T> &s) { return s.data(); }
template <typename T>
T *__getPtr(gsl::span<T> &s) { return s.data(); }
/*
 * Access class used by BaseConnectivity functions in order to free Connectivity from its link to element set
*/
template<typename A>
class DirectAccess {
public :
	static int getNum(const A* oc, int i)
	{return oc->num(i); }
	static const int * getData(const A* oc, int i)
	{return __getPtr( (*oc)[i] );}
	static int getNumTarget(const A *oc)
	{return oc->getNumTarget();}
	/*static int * getTarget(A* oc)
	  {return oc->getTarget();}*/
	/*static int * getPointer(A* oc)
	  {return oc->getPointer();}*/
	static int getSize(const A* oc)
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
	Connectivity reverse() const;
	Connectivity* altReverse(float *w = nullptr);
	template <class B, class AB>
	Connectivity* transcon(const BaseConnectivity<B, AB> *tc) const;
	template <class B, class AB>
	Connectivity transcon(const BaseConnectivity<B,AB>& tc) const;

	template <class B, class AB>
	Connectivity* altTranscon(const BaseConnectivity<B,AB>& tc) const;
	/// \brief Transitive connectivity but only where owned[ tc[i][j] ] == me.
	template <class B, class AB>
	Connectivity* transcon(const BaseConnectivity<B, AB> *tc, int *ownedBy, int me);

	/* following functions are used to make transcon accept a BaseConnectivity as argument */
//    int *operator[](int i) {return(Accessor::getData(static_cast<A*>(this),i));}
//	const int *operator[](int i) const {return(Accessor::getData(static_cast<const A*>(this),i));}
    gsl::span<const int> operator[](int i) const {
    	return { Accessor::getData(static_cast<const A*>(this),i), Accessor::getNum(static_cast<const A*>(this), i) };
    }
	int csize() const {return(Accessor::getSize(static_cast<const A*>(this)));}
	int num(int i) const {return(Accessor::getNum(static_cast<const A*>(this),i));}
	int getNumTarget() const {return(Accessor::getNumTarget(static_cast<const A*>(this)));}
	void print() const
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

template <typename T>
int getNumNodes(const T &t) { return t.numNodes(); }
template <typename T>
auto getNodes(const T &t, int *nd) { return t.nodes(nd); }

/** Class to access a set of element's individual connections */
template <class A>
class SetAccess {
public:
	SetAccess(const A &set) : sz(set.last()), set(set) {}
	SetAccess(const A &set, int sz) : sz(sz), set(set) {}
	int size() const { return sz; }; //<! returns the number of members of the set
	int numNodes(int i) const { return (set[i] ? getNumNodes(*set[i]) : 0); } //<! returns the number of targets for member i
	void nodes(int i, int *nd) const { if(set[i]) getNodes(*set[i], nd); }; //<! copies into nd the targets for member i
private:
	int sz;
	const A &set;
};

class Connectivity : public BaseConnectivity<Connectivity,DirectAccess<Connectivity> >
{
protected:
	int size;           // size of pointer
	std::vector<int> pointer;       // pointer to target
	std::vector<int> target;        // value of the connectivity
	std::vector<float> weight;      // weights of pointer (or 0)

public:
	typedef int IndexType;
	int getNumTarget() const {return target.size(); }
	int * getTarget() {return target.data(); }
	const int * getTarget() const {return target.data(); }
        int * getPointer() {return pointer.data(); }
	const int * getPointer() const {return pointer.data(); }

	/** \brief Factory static method.
	 *
	 * \details The targets are kept in the order in which they are given in the data.
	 * Such property is important for building consistent communication lists.
	 *
	 * @tparam RangeT A type on which one can iterate to obtain
	 * pair like objects of source, target named first, second.
	 * @param range Range object.
	 * @return
	 */
	template <typename RangeT>
	static Connectivity fromLinkRange(const RangeT &range);


	Connectivity() { size = 0; }
	/** \brief Constructor for any object that is equipped with the methods of a set.
	 *
	 * @tparam A A type sastisfying the notion of a set.
	 * @param sa The set.
	 */
	template <class A>
	Connectivity(const SetAccess<A> &sa);
	/** \brief Build a connectivity which is a map from SommerElement to the actual (single) element to which
	 * the sommerElement is attached. TODO Get rid of this. The result is a one to one map. Not a connectivity.
	 */
	Connectivity(Elemset *, int, SommerElement**, bool wetFlag = false);
	Connectivity(int _size, int *_pointer, int *_target, int _removeable=1, float *_weight = 0);
	Connectivity(int _size, std::vector<int> _pointer, std::vector<int> _target,
		             std::vector<float> _weight = std::vector<float>{});
	Connectivity(int _size, int *count);
	Connectivity(int _size, int count);
	Connectivity(BinFileHandler &, bool oldSower = false);
	Connectivity(int ns); //dec
	Connectivity(FILE *f, int nElem); // JAT 100614
	/** \brief Construct a connectivity accounting for Lagrange multipliers.
	 *
	 * @param els
	 * @param nodeToElem
	 */
	Connectivity(const Elemset &els, Connectivity *nodeToElem);
	Connectivity(const Connectivity&) = default;
	Connectivity(Connectivity &&) = default;
	Connectivity &operator=(Connectivity &&) = default;
	Connectivity &operator=(const Connectivity &) = default;
	size_t write(BinFileHandler& f);
	size_t writeg(BinFileHandler& f);
	size_t read(FILE* f);
        size_t write(FILE* f);

	void countlink(int from, int to); //DEC
	void addlink(int from, int to); //DEC

//	Connectivity(SommerElement  **, int);
	Connectivity(SommerElement  **, int, int);

	virtual ~Connectivity();
	virtual void end_count(); //dec
	gsl::span<const int> operator[](int i) const;
	gsl::span<int> operator[](int i);
	/// \brief Obtain the number of source nodes.
	int csize() const;
	int numConnect() const; // Total number of connections
	/// \brief Get all targets.
	auto &allTargets() const { return target; }
	int num(int) const;
	int num(int, int*) const;
	int offset(int i)  const { return pointer[i]; } // Begining of ith part
	int offset(int i,int j) const; // returns a unique id for connection i to j
	int cOffset(int i,int j) const; // returns offset(i,j) - pointer[i]
	bool locate(int i,int j) const
	{ for(int k=pointer[i]; k<pointer[i+1]; ++k) if(target[k] == j) return true; return false; }

	Connectivity* alloc_reverse() const
	{ return new Connectivity{BaseConnectivity<Connectivity>::reverse()}; } // creates t->s from s->t

	Connectivity* transconOne(Connectivity*);
	int getTargetValue(int i) { return target[i]; }

	void findPseudoDiam(int *n1, int *n2, int *mask=0);
	int  rootLS(int root, int *xls, int *ls, int &w, int *mask=0);

	long long memsize() {return ((long long)size + getNumTarget() + 1)*sizeof(int);}

	// Create a rooted level structure
	int *renumSloan(int *mask, int &firstNum, int *ren = 0);
	int *renumRCM(int *mask, int &firstNum, int *ren = 0);
	compStruct renumByComponent(int);
	void print(FILE * = stderr, int node=-1) const;
	int findMaxDist(int *) const;
	int findProfileSize(EqNumberer *eqNumber, int unroll=1) const;
	const std::vector<int> &ptr() const { return pointer; }
	auto &tgt() { return target; }
	auto &tgt() const { return target; }
	/** \brief Create a connectivity that appends a given connectivity at the end of this one.
	 * \details The number of sources is the sum of both. */
	Connectivity *alloc_append(Connectivity *cn) const;
	// Create a copy of this connectivity without...
	Connectivity *collapse();
	Connectivity *subSection(bool *);
	Connectivity *trim(Connectivity *);
	void sortTargets();
	void renumberTargets(int *map);
	void renumberTargets(std::map<int, int> &);
	int numNonZeroP() const;

	bool isDiagonal() const; // returns true if each target is only connected to itself
	Connectivity modify();
	/// \brief Create a similar connectivity with a self connection for nodes that did not have one.
	Connectivity withSelfConnection() const;
	Connectivity *modifyAlt();
	void combine(const Connectivity *con2, std::vector<int> &cmap, const std::vector<int> &cmap2);  // adds con2 to this
	// adds all the entries in cmap (of size addSize)to each of the line in the current connectivity specified by entries in cmap
	// e.g. cmap = [2 3] and (*this)[2] = [1 2 4 5] (*this)[3] = [3 5], then (*this)[2] becomes [1 2 4 5 3]; (*this)[3] becomes [3 5 2]
	Connectivity combineAll(int addSize, int *cmap);

	double estimateComponentCost(EqNumberer *eqn, compStruct &cs, double *cost, double *bandwidth,
	                             double coef=400, int unroll=1);
	double estimateCost(EqNumberer *eqn, double &cost, double &bandwidth,
	                    double coef=400, int unroll=1);

	Connectivity * SCOTCH_graphPart(int partnbr);
};



inline int
Connectivity::csize() const { return size; }

inline int
Connectivity::num(int n) const { return (n < size) ? pointer[n+1] - pointer[n] : 0; }

inline gsl::span<const int>
Connectivity::operator[](int i) const { return { target.data() +pointer[i], num(i) }; }

inline gsl::span<int>
Connectivity::operator[](int i) { return { target.data() +pointer[i], num(i) }; }

inline int
Connectivity::numConnect() const { return getNumTarget(); }

inline bool
Connectivity::isDiagonal() const { return (getNumTarget() == size) ? true : false; }

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
	mutable int cacheVal, cacheSize;
	mutable int *cacheTg;
	mutable int cacheNumVal, cacheNumIdx;
public:
	ImplicitConnectivity(A _a) : a(_a) { cacheSize = 0; cacheTg = 0;
		cacheVal=-1; cacheNumIdx = -1; }
	~ImplicitConnectivity() { if(cacheTg) delete [] cacheTg; }
	int csize() const { return Accessor::getSize(*a); }
	int num(int i) const { if(cacheNumIdx != i) {
			cacheNumVal = Accessor::getNum(*a,i);
			cacheNumIdx = i;
		}
		return cacheNumVal;
	}
	int *getNewCacheVal(int j) const {
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
	const int *operator[](int i) const {
		if(cacheVal == i)
			return cacheTg;
		else
			return getNewCacheVal(i);
	}
	int getNumTarget() const {
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
Connectivity
BaseConnectivity<A,Accessor>::reverse() const
{
	// The reverse connectivity has the same size as the original
	int size = csize(); //Accessor::getSize(static_cast<A*>(this));
	int numTarget = getNumTarget(); //Accessor::getNumTarget(static_cast<A*>(this));
	std::vector<int> res_target(numTarget);

	// Find the max of target
	int maxtarg = -1; // PJSA
	int i;
	for(i=0; i < size; ++i) {
		int nTg = Accessor::getNum(static_cast<const A*>(this), i);
		auto tg = Accessor::getData(static_cast<const A*>(this), i);
		for(int j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
	}
	int res_size = maxtarg+1;
	std::vector<int> res_pointer(res_size+1);
	for(i = 0; i <= res_size; ++i) res_pointer[i] = 0;

	// Now do a first pass to fill in res_pointer
	for(i=0; i < size; ++i) {
		int nTg = Accessor::getNum(static_cast<const A*>(this), i);
		auto tg = Accessor::getData(static_cast<const A*>(this), i);
		for(int j = 0; j < nTg; ++j) res_pointer[tg[j]]++;
	}

	for(i = 1; i <= res_size; ++i) res_pointer[i] += res_pointer[i-1];

	// Second pass fills in target
	for(i=0; i < size; ++i) {
		int nTg = Accessor::getNum(static_cast<const A*>(this), i);
		auto tg = Accessor::getData(static_cast<const A*>(this), i);
		for(int j = 0; j < nTg; ++j)
			res_target[--res_pointer[tg[j]]] = i;
	}


	return Connectivity(res_size, std::move(res_pointer), std::move(res_target));
}

// Important NOTE: transcon cannot be called with the tc == this if tc is
// an Implicit Connectivity!!!
template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::transcon(const BaseConnectivity<B, AB> *tc) const
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
	std::vector<int>flags(tgmax, -1);

	// Compute the new pointers
	std::vector<int> np(size+1);
	int cp = 0;
	for(i = 0; i < size; ++i) {
		np[i] = cp;
		int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
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
	flags.assign(tgmax, -1);
	std::vector<int> ntg(cp);
	cp = 0;
	for(i = 0; i < size; ++i) {
		int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
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
	return new Connectivity(size, std::move(np), std::move(ntg));
}

// Important NOTE: transcon cannot be called with the tc == this if tc is
// an Implicit Connectivity!!!
template<typename A, class Accessor>
template<class B, class AB>
Connectivity BaseConnectivity<A,Accessor>::transcon(const BaseConnectivity<B,AB>& tc) const
{
	int i,j,k;

	// First find the biggest target so we can size arrays correctly
	int tgmax=-1;

	int size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
	for(i = 0; i < tc.csize(); ++i)
		for(j = 0; j < tc.num(i); ++j)
			tgmax = std::max(tgmax, tc[i][j]);
	tgmax++; // Important adjustment

	// Now we can size the array that flags if a target has been visited
	std::vector<int>flags(tgmax, -1);

	// Compute the new pointers
	std::vector<int> np(size+1);
	int cp = 0;
	for(i = 0; i < size; ++i) {
		np[i] = cp;
		int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
		for(j = 0; j < nTg; ++j) {
			int intermed = tg[j];
			for(k = 0; k < tc.num(intermed); ++k)
				if(flags[tc[intermed][k]] != i) {
					flags[tc[intermed][k]] = i;
					cp ++;
				}
		}
	}
	np[size] = cp;

	// Now allocate and fill the new target
	flags.assign(tgmax, -1);
	std::vector<int> ntg(cp);
	cp = 0;
	for(i = 0; i < size; ++i) {
		int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
		for(j = 0; j < nTg; ++j) {
			int intermed = tg[j];
			for (k = 0; k < tc.num(intermed); ++k)
				if(flags[tc[intermed][k]] != i) {
					flags[tc[intermed][k]] = i;
					ntg[cp] = tc[intermed][k];
					cp ++;
				}
		}
	}
	return { size, std::move(np), std::move(ntg) };
}

template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::transcon(const BaseConnectivity<B, AB> *tc, int *ownedby, int me)
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
	std::vector<int>flags(tgmax, -1);

	// Compute the new pointers
	int *np = new int[size+1];
	int cp = 0;
	for(i = 0; i < size; ++i) {
		np[i] = cp;
		int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
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
	flags.assign(tgmax, -1);
	int *ntg = new int[cp];
	cp = 0;
	for(i = 0; i < size; ++i) {
		int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
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
		auto tg = Accessor::getData(static_cast<A*>(this), i);
		for(int j = 0; j < nTg; ++j) maxtarg = std::max(tg[j],maxtarg);
	}
	int res_size = maxtarg+1;
	int *res_pointer = new int[res_size+1];
	for(i = 0; i <= res_size; ++i) res_pointer[i] = 0;

	// Now do a first pass to fill in res_pointer
	for(i=0; i < size; ++i) {
		int nTg = Accessor::getNum(static_cast<A*>(this), i);
		auto tg = Accessor::getData(static_cast<A*>(this), i);
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
		auto tg = Accessor::getData(static_cast<A*>(this), i);
		for(int j = 0; j < nTg; ++j)
			res_target[res_pointer[tg[j]] + count[tg[j]]++] = i;
	}
	delete [] count;

	Connectivity *res = new Connectivity(res_size, res_pointer, res_target, 1, w);
	return res;
}

template<typename A, class Accessor>
template<class B, class AB>
Connectivity* BaseConnectivity<A,Accessor>::altTranscon(const BaseConnectivity<B,AB>& tc) const
{
	// PJSA 12-12-05 this version of transcon doesn't include connections with self unless entirely internal
	int i,j,k;

	// First find the biggest target so we can size arrays correctly
	int tgmax=-1;

	int size = csize(); //Accessor<A>::getSize(static_cast<A*>(this));
	for(i = 0; i < tc.csize(); ++i)
		for(j = 0; j < tc.num(i); ++j)
			tgmax = std::max(tgmax, tc[i][j]);
	tgmax++; // Important adjustment

	// Now we can size the array that flags if a target has been visited
	std::vector<int> flags(tgmax);
	for(i = 0; i < tgmax; ++i) flags[i] = -1;

	// Compute the new pointers
	std::vector<int> np(size+1);
	int cp = 0;
	for(i = 0; i < size; ++i) {
		np[i] = cp;
		int nTg =  num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
		for(j = 0; j < nTg; ++j) {
			int intermed = tg[j];
			for(k = 0; k < tc.num(intermed); ++k) {
				if((tc[intermed][k] == i) && (tc.num(intermed) != 1)) continue;
				if(flags[tc[intermed][k]] != i) {
					flags[tc[intermed][k]] = i;
					cp ++;
				}
			}
		}
	}
	np[size] = cp;

	// Now allocate and fill the new target
	for(i = 0; i < tgmax; ++i)
		flags[i] = -1;
	std::vector<int> ntg(cp);
	cp = 0;
	for(i = 0; i < size; ++i) {
		int nTg = num(i); //Accessor<A>::getNum(static_cast<A*>(this), i);
		auto tg = (*this)[i]; //Accessor<A>::getData(static_cast<A*>(this), i);
		for(j = 0; j < nTg; ++j) {
			int intermed = tg[j];
			for (k = 0; k < tc.num(intermed); ++k) {
				if((tc[intermed][k] == i) && (tc.num(intermed) != 1)) continue;
				if(flags[tc[intermed][k]] != i) {
					flags[tc[intermed][k]] = i;
					ntg[cp] = tc[intermed][k];
					cp ++;
				}
			}
		}
	}

	Connectivity *res = new Connectivity(size, np, ntg);
	return res;
}

template <class A>
Connectivity::Connectivity(const SetAccess<A> &sa)
{
	int i;

	size = sa.size();

	// Find out the number of targets we will have
	pointer.resize(size+1) ;
	int pp = 0;
	for(i=0; i < size; ++i) {
		pointer[i] = pp;
		pp += sa.numNodes(i);
	}
	pointer[size] = pp;

	// Create the target array
	target.resize(pp);

	// Fill it in
	for(i=0; i < size; ++i) {
		sa.nodes(i, target.data()+pointer[i]);
	}
}

template<typename RangeT>
Connectivity Connectivity::fromLinkRange(const RangeT &range) {
	auto map = [] (int n) { return n; };
	decltype(map(range.begin()->first)) maxIdx{0};
	for(auto &p : range) {
		auto idx = map(p.first);
		maxIdx = std::max(maxIdx, idx);
	}
	auto size = maxIdx+1;
	std::vector<int> pointers(size+1, 0);
	for(auto &p : range)
		++pointers[map(p.first)];
	for(size_t i = 0; i < size; ++i)
		pointers[i+1]+=pointers[i];
	std::vector<int> targets;
	targets.resize(pointers[size]);
	for(auto &p : range)
		targets[--pointers[map(p.first)]] = p.second;
	return Connectivity(size, std::move(pointers), std::move(targets));
}
#endif
