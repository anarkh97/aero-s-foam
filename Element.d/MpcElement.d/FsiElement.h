#ifndef _FSIELEMENT_H_
#define _FSIELEMENT_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
class LMPCons;

class FsiElement : public Element
{
	int nnodes;  // number of nodes
	int *nn;  // node numbers
	int ndofs; // number of dofs
	int *renumTable;
	LMPCons *fsi;

public:
	explicit FsiElement(LMPCons *mpc);

	~FsiElement() override { delete [] nn; delete [] renumTable; };

	void setProp(StructProp *p) { };
	bool isFsiElement() override { return true; }
// JLchange: only valid for the case of one-one connection fsi
	int fsiFluidNode() override { return nn[nnodes-1]; }
	int fsiStrutNode() override { return nn[0]; }

	LMPCons* cons() { return fsi; }

	void renum(int *table) override;
	void renum(EleRenumMap&) override;

	FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg=1) const override;
	FullSquareMatrix imagStiffness(CoordSet&, double *kel, int flg=1);
	FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;

	void markDofs(DofSetArray &) override;
	int* dofs(DofSetArray &, int *p) override;
	int  numDofs() const override { return ndofs; }
	int  numNodes() const override { return nnodes; }
	int* nodes(int *) const override;
	int  numInternalNodes() override { return 0; }

	int  getTopNumber() override { return 502; }
	int  numTopNodes() override { return nnodes; }

	PrioInfo examine(int sub, MultiFront *mf) override;
// JLchange    bool isSafe() const override { return false; }
};
#endif

