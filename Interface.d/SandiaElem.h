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
  FullSquareMatrix stiffness(const CoordSet&, double *, int = 1) const override;
  FullSquareMatrix massMatrix(const CoordSet&, double *, int = 1) const override;
  void markDofs(DofSetArray &dsa);
  int* dofs(DofSetArray &dsa, int *p=0);
   int numDofs() const override;
  int numNodes() const override;
  int* nodes(int * = 0) const override;
  bool isSafe() const;
  bool isRotMidSideNode(int iNode);
};

