#include <Element.d/Element.h>

class SandiaElem : public Element
{
  int nnodes, *nds;
  unsigned short int *flags;
  Element *proxy;

 public:
  SandiaElem(int numberOfNodes, int* nodeNumbers, unsigned short int* dofFlags);
  ~SandiaElem() { if(nds) delete [] nds; }
  void renum(int *renumberingTable);
  FullSquareMatrix stiffness(CoordSet&, double *, int = 1);
  FullSquareMatrix massMatrix(CoordSet&, double *, int = 1);
  void markDofs(DofSetArray &dsa);
  int* dofs(DofSetArray &dsa, int *p=0);
  int numDofs();
  int numNodes();
  int* nodes(int * = 0);
  bool isSafe();
  bool isRotMidSideNode(int iNode);
};

