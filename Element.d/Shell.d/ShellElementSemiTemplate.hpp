#ifndef _SHELLELEMENTSEMITEMPLATE_HPP_
#define _SHELLELEMENTSEMITEMPLATE_HPP_

#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

template<typename doublereal>
class ShellElementSemiTemplate 
{
  public:
    void sands8(doublereal *_xl,doublereal *_yl, doublereal *_zl,
         doublereal e, doublereal nu, doublereal *_h,
         doublereal *_v, doublereal *_stress,
         bool strainFlg, int surface, doublereal thrmStr);
    
    void tria3d(bool flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
                doublereal e, doublereal nu, doublereal *_h, doublereal *_rk);

    void tria3dthickness(bool flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
                         doublereal e, doublereal nu, doublereal *_h, doublereal *_drkdthick);

    void vms8WRTdisp(doublereal *_xl,doublereal *_yl, doublereal *_zl,
                     doublereal e, doublereal nu, doublereal *_h,
                     doublereal *_v, doublereal *_vmsWRTdisp,
                     bool strainFlg, int surface, doublereal thrmStr = 0);

    void vms8WRTthic(doublereal *_xl,doublereal *_yl, doublereal *_zl,
                     doublereal e, doublereal nu, doublereal *_h,
                     doublereal *_v, doublereal *_vmsWRthic,
                     bool strainFlg, int surface, doublereal thrmStr = 0);

  protected:
    virtual void rotation(doublereal *_v1o, doublereal *_v2o, doublereal *_v3o,
                          doublereal *_v1n, doublereal *_v2n, doublereal *_v3n,
                          doublereal *_r);

    virtual void membra(doublereal *_x, doublereal *_y, doublereal alpha,
                        int *_le, doublereal *_q);

    virtual void transform(doublereal *_xl, doublereal *_yl, doublereal *_zl,
                           doublereal *_xg, doublereal *_yg, doublereal *_zg,
                           doublereal *_str, doublereal *_tMat);

    virtual void vonmis(doublereal rmx, doublereal rmy, doublereal rmxy,
                        doublereal rnx, doublereal rny, doublereal rnxy,
                        doublereal &t,  doublereal &sv, int surface=1);

    virtual void compj2(doublereal sx, doublereal sy, doublereal sxy, doublereal &svm);

    void momen(doublereal *_x, doublereal *_y, int *_lb, doublereal *_l);

    void straineq(doublereal rmx, doublereal rmy, doublereal rmxy,
                  doublereal rnx, doublereal rny, doublereal rnxy,
                  doublereal nu,  int surface, doublereal &ebar);

    void equiv(doublereal ex, doublereal ey, doublereal ez, doublereal exy, doublereal &eq);

    void basico(doublereal *_x, doublereal *_y, doublereal *_db, doublereal f,
                doublereal clr, doublereal cqr, int *_ls, doublereal *_sm, int m, char *status);

    void sm3mb(doublereal *_x, doublereal *_y, doublereal *_dm,
               doublereal alpha, doublereal f, int *_ls, doublereal *_sm, int m, char *status);

    void smcbh(doublereal *_x, doublereal *_y, doublereal *_db,
               doublereal f, int *_ls, doublereal *_sm, int m, char *status);
  
    void sm3mhe(doublereal *_x, doublereal *_y, doublereal *_dm, doublereal f, int *_ls, doublereal *_sm, int m, char *status);

    void trirotation(doublereal *_sm, doublereal *_r1); 
};

#endif
