#ifndef _SHELLELEMENTSEMITEMPLATE_HPP_
#define _SHELLELEMENTSEMITEMPLATE_HPP_

#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

// Note: the following needs to be defined for full compatability with the implementation of elements 8 and 20
//#define COMPATABILITY_MODE

template<typename doublereal>
class ShellElementSemiTemplate 
{
  public:
    void sands8(doublereal *_xl,doublereal *_yl, doublereal *_zl,
         doublereal e, doublereal nu, doublereal *_h,
         doublereal *_v, doublereal *_stress,
         bool strainFlg, int surface, doublereal thrmStr);

    void vms8WRTdisp(doublereal *_xl,doublereal *_yl, doublereal *_zl,
                     doublereal e, doublereal nu, doublereal *_h,
                     doublereal *_v, doublereal *_vmsWRTdisp,
                     bool strainFlg, int surface, doublereal thrmStr = 0);

    void vms8WRTthic(doublereal *_xl,doublereal *_yl, doublereal *_zl,
                     doublereal e, doublereal nu, doublereal *_h,
                     doublereal *_v, doublereal *_vmsWRthic,
                     bool strainFlg, int surface, doublereal thrmStr = 0);

  protected:
    void rotation(doublereal *_v1o, doublereal *_v2o, doublereal *_v3o,
                  doublereal *_v1n, doublereal *_v2n, doublereal *_v3n,
                  doublereal *_r);

    void membra(doublereal *_x, doublereal *_y, doublereal alpha,
                int *_le, doublereal *_q);

    void transform(doublereal *_xl, doublereal *_yl, doublereal *_zl,
                   doublereal *_xg, doublereal *_yg, doublereal *_zg,
                   doublereal *_str, doublereal *_tMat);

    void vonmis(doublereal rmx, doublereal rmy, doublereal rmxy,
                doublereal rnx, doublereal rny, doublereal rnxy,
                doublereal &t,  doublereal &sv, int surface);

    void compj2(doublereal sx, doublereal sy, doublereal sxy, doublereal &svm);

    void momen(doublereal *_x, doublereal *_y, int *_lb, doublereal *_l);

    void straineq(doublereal rmx, doublereal rmy, doublereal rmxy,
                  doublereal rnx, doublereal rny, doublereal rnxy,
                  doublereal nu,  int surface, doublereal &ebar);

    void equiv(doublereal ex, doublereal ey, doublereal ez, doublereal exy, doublereal &eq);

};

#endif
