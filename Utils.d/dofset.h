#ifndef _DOFSET_H_
#define _DOFSET_H_
#include <iostream>

class DofSet {
  protected:
    int flags;
  public:
    static const int max_known_dof = 24;
    static const int max_known_nonL_dof = 9;
    static int 
      Xdisp,
      Ydisp,
      Zdisp,
      XYZdisp,
      Xrot,
      Yrot,
      Zrot,
      XYZrot,
      Temp,
      Helm,
      Contact,
      IntPress,
      Potential,
      LagrangeE,
      LagrangeI;
    static int nonL_dof;

    static const int DispAndRot = 0x3f;
    static DofSet nullDofset;

    // constructors
    DofSet()      { flags = 0; }
    DofSet(int t) { flags = t; }

    DofSet & operator |= (const DofSet &ds) { flags |= ds.flags; return *this; }
    bool operator == (const DofSet &ds) const { return flags == ds.flags; }

    // mark marks given dofs as being used
    void mark(int dofs) { flags |= dofs; }

    void unmark(int dofs) { //fprintf(stderr,"In unmark: ");
                            //fprintf(stderr,"flags = %d; ", flags);
                            //fprintf(stderr,"dofs = %d\n", dofs);
                            flags = flags ^ (flags & dofs); }
    void unmarkAll() { flags = 0; }

    // tests if all the given dofs are present
    int test(int dofs) { return (flags & dofs) == dofs ; }

    // see presence of at least one dof
    bool contains(int dofs) { //fprintf(stderr,"In contains: ");
                              //fprintf(stderr,"flags = %d; ", flags);
                              //fprintf(stderr,"dofs = %d\n", dofs);
                              return (flags & dofs) != 0; }

/*
    bool containsAllDisp(int dim) {
      if(dim >=3) return contains(XYZdisp);
      else if(dim == 2) return contains(Xdisp & Ydisp);
      else if(dim == 1) return contains(Xdisp);
      else return true;
    }
*/
    bool containsAllDisp(int dim) {
      switch(dim) {
        case 1 :
          return contains(Xdisp);
          break;
        case 2:
          return (contains(Xdisp) && contains(Ydisp));
          break;
        case 3:
          return (contains(Xdisp) && contains(Ydisp) && contains(Zdisp));
          break;
        default:
          std::cerr << " *** WARNING: in DofSet::containsAllDisp(), dim = " << dim << std::endl;
          return false;
          break;
      }
    }

    bool containsAnyRot() {
      return contains(Xrot | Yrot | Zrot);
    }

    // counts the number dofs marked in this node
    int count();

    int number(DofSet, int *); // This routine complements the previous
                  // one, numbering all the DOFs in the first argument

    // Locates a dof in the current set. Only works for single dofs
    // i.e. result is undefined for composites such as XYZrot and XYZdisp
    int locate(int dof) const;

    int list() { return flags; }

    DofSet operator & (DofSet d) { return DofSet(flags & d.flags); }
    DofSet operator | (DofSet d) { return DofSet(flags | d.flags); }
    DofSet operator ^ (DofSet d) { return DofSet(flags ^ d.flags); }
    DofSet & operator &=(const DofSet &d) { flags &= d.flags; return *this; } //HB
    void print(char* msg=0); //HB
};

class Elemset;
class Element;

class EqNumberer {
  protected:
    int numnodes; 	// total number of nodes
    int *node_offset; 	// Where the nodes are
    int *node_num_dofs; // Number of dofs associated with a node
    int *renummap;    	// Renumbering mapping of the nodes
    int myMap;

  public:
    EqNumberer() { node_offset = 0; node_num_dofs = 0; renummap = 0; myMap = 0; }
    virtual ~EqNumberer();

    // Return the total number of degrees of freedom
    int size() { return (node_offset) ? node_offset[numnodes] : 0; }

    // Return the number of nodes
    int numNodes() const { return numnodes; }

    // Return the first dof of a node
    int firstdof(int node) 
     { return (node_num_dofs[node] > 0) ? node_offset[node] : -1; }

    // Return the weight of a node (its number of dofs)
    int weight(int n) { return node_num_dofs[n]; }

    int *allWeights() { return node_num_dofs; }

    int *allOffsets() { return node_offset; }
    int *renumPtr() { return renummap; }
    void print();
};

class DofSetArray : public EqNumberer {
  protected:
    DofSet *dofs;
    int *rowcolnum;
    int *invrowcol;
    int *dofType; // 0 = translational, 1 = rotational
    bool myDofs; 
    DofSetArray() { dofs = 0; rowcolnum = 0; invrowcol = 0; dofType = 0; myDofs = true; }

  protected:
    void  makeOffset();
    void  makeModifiedOffset();
  public:
    DofSetArray(int nnode, int *dofsPerNode, int *renumtable); // for DEC
    DofSetArray(int nnodes, Elemset &elearray, int *renumtable=0, int myMap=0);
    DofSetArray(int nnodes, int *renumtable=0, int myMap=0);
    DofSetArray(Element *ele);

    virtual ~DofSetArray();

    void initialize() { dofs = 0; rowcolnum = 0; invrowcol = 0; dofType = 0; myDofs = false; }

    // locate a dof for a given node
    int locate(int node, int dof) const;
    int number(int node, DofSet, int *);

    // Mark dofs for a node
    void mark(int node, int dof);
    void mark(int *node, int numNode, int dof);

    // Return the DofSet of a node
    DofSet &operator [](int i) { return dofs[i]; }

    int getRCN(int dof)   { return rowcolnum[dof]; }
    int invRCN(int dof)   { return invrowcol[dof]; }

    int* getUnconstrNum() { return rowcolnum; }
    int* getConstrndNum() { return invrowcol; }

    int* makeDofTypeArray(); // creates and returns dofType array
                             // 0 = translational, 1 = rotational
    int* getDofTypeArray() { return dofType; }

    void setWeight(int n, int w); 
    int getWeight(int n);
    void finish() { makeModifiedOffset(); }
    friend class ConstrainedDSA;
    void clean_up();
};

template <typename VecType>
void
zeroRotDofs(const DofSetArray &dsa, VecType &vec) {
  static const int ROT_DOFS[] = { DofSet::Xrot, DofSet::Yrot, DofSet::Zrot };
  static const int ROT_DOFS_SIZE = sizeof(ROT_DOFS) / sizeof(ROT_DOFS[0]);

  const int nodeCount = dsa.numNodes();
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    for (const int *dofType = ROT_DOFS; dofType != ROT_DOFS + ROT_DOFS_SIZE; ++dofType) {
      const int dofLoc = dsa.locate(iNode, *dofType);
      if (dofLoc >= 0) {
        vec[dofLoc] = 0.0;
      }
    }
  }
}

class BCond;
class ComplexBCond;

class ConstrainedDSA : public DofSetArray {
  public:
    ConstrainedDSA(DofSetArray &, int, BCond *);
    ConstrainedDSA(DofSetArray &, DofSetArray &, int, BCond *, int);
    ConstrainedDSA(DofSetArray &, int, BCond *, int *bc);
    ConstrainedDSA(DofSetArray &, int, ComplexBCond *, int *bc);
    ConstrainedDSA(DofSetArray &dsa, BCond *bcdata, int nbc,
                   ComplexBCond *bcd, int nbcd, int *bc);

    ConstrainedDSA(DofSetArray &dsa, int nbc, BCond *bcond,
                   int numCornerNodes, int *cornerNodes, DofSet *cornerDofs,
                   int ncbc = 0, ComplexBCond *cbcond = 0, int numWetInterfaceNodes = 0,
                   int *wetInterfaceNodes = 0, DofSet *wetInterfaceDofs = 0);

    ConstrainedDSA(DofSetArray &dsa, int n, int *sing = 0);  // PJSA: 1-23-01
    ConstrainedDSA(DofSetArray &dsa, ConstrainedDSA &cdsa, int n, int *sing = 0);
    ConstrainedDSA(DofSetArray &dsa, ConstrainedDSA &cdsa);
    virtual ~ConstrainedDSA() { /* nothing to delete */ };
    int getInvRCNmax() { return invrowcolmax; }  // PJSA: 11-12-02
  private:
    int invrowcolmax;
};

// Auxiliary definitions
#define BCFREE  0
#define BCLOAD  1
#define BCFIXED 2

class SimpleNumberer : public EqNumberer {
   public:
      SimpleNumberer(int nnodes, int *renumb = 0, int myMap=0);
      virtual ~SimpleNumberer() { /* nothing to delete */ };
      void setWeight(int, int);
      void makeOffset();
};

#endif
