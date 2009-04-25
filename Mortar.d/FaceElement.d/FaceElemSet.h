// ---------------------------------------------------------------- 
// HB - 06/25/03
// ---------------------------------------------------------------- 
// WARNING: IT IS CURRENTLY IMPLICITLY ASSUMED THAT NO GAP EXIST IN
//          THE NUMBERING OF THE FACE ELEMENTS IN THE FaceElemSet,
//          AND THE NUMBERING START FROM 0 TO N-1 (WHERE N IS THE 
//          NUMBER OF FACE ELEMENTS REALLY INSTANCIATED), SO THAT
//          FaceElemSet::nElems() RETURN N.
// ---------------------------------------------------------------- 
#ifndef _FACEELEMSET_H_
#define _FACEELEMSET_H_

// STL
#include <map>

// FEM headers
#include <Utils.d/BlockAlloc.h>
//#include <Mortar.d/FaceElement.d/FaceElement.h>

#ifdef SOWER_SURFS
#include <Utils.d/BinFileHandler.h>
#endif
                                                                                                              
class FaceElement;

class FaceElemSet {
  protected:
    FaceElement **elem;
    int emax;
    int nPhantoms;
    BlockAlloc ba;
  public:
    FaceElemSet(int = 256);
    //~FaceElemSet() { ba.~BlockAlloc(); deleteElems(); }
    ~FaceElemSet() { deleteElems(); }
    int size() { return emax; }
    int last();
    int numPhantoms()  { return nPhantoms; }
    FaceElement *operator[] (int i) { return elem[i]; }
    void elemadd(int num, FaceElement *);
    void elemadd(int num, int type, int nnodes, int *nodes);
    void setEmax(int max)  { emax = max; }
    void setNumPhantoms(int n)  { nPhantoms = n; }
    void list();

    int nElems();
    void print();

    void Renumber(std::map<int,int>& OldToNewNodeIds);

    void deleteElems()  { if(elem){ delete [] elem; elem = 0; } emax = 0; nPhantoms = 0; }

#ifdef SOWER_SURFS
    void WriteSower(BinFileHandler& file, int* ndMap=0);
#endif
};

#endif
