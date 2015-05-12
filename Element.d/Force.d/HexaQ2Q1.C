#if defined(USE_EIGEN3) && (__cplusplus >= 201103L) && defined(HAS_CXX11_TEMPLATE_ALIAS)
#include <Element.d/Force.d/HexaQ2Q1.h>

const DofSet HexaQ2Q1::NODALDOFS[20] = { DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp,
                                         DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp, DofSet::XYZdisp };

PrioInfo examineHex20(int sub, MultiFront *mf, int *nn);

HexaQ2Q1::HexaQ2Q1(int* _nn)
 : MixedFiniteElement<HexaQ2Q1ThreeFieldStrainEnergyFunction>(20, const_cast<DofSet*>(NODALDOFS), _nn)
{}

int
HexaQ2Q1::getTopNumber()
{
  return 172;
}

PrioInfo
HexaQ2Q1::examine(int sub, MultiFront *mf)
{
  return examineHex20(sub, mf, nn);
}

int
HexaQ2Q1::getQuadratureOrder()
{
  return 3;
}

#endif
