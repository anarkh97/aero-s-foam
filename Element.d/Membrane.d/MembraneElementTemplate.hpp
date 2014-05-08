#ifndef _MEMBRANEELEMENTTEMPLATE_HPP_
#define _MEMBRANEELEMENTTEMPLATE_HPP_

#include <Element.d/Shell.d/ShellElementSemiTemplate.hpp>

template<typename doublereal>
class MembraneElementTemplate: public ShellElementSemiTemplate<doublereal> 
{
  public:
     MembraneElementTemplate() {}
    ~MembraneElementTemplate() {}
  
     void sands19(doublereal *_xl, doublereal *_yl, doublereal *_zl,
                  doublereal e, doublereal nu, doublereal *_h,
                  doublereal *_v, doublereal *_stress,
                  bool strainFlg);
     void trimem(int flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
                 doublereal e, doublereal nu, doublereal *_h, doublereal *_rk); 
/*
  protected:
     void rotation(doublereal *_v1o, doublereal *_v2o, doublereal *_v3o,
                   doublereal *_v1n, doublereal *_v2n, doublereal *_v3n,
                   doublereal *_r);

     void membra(doublereal *_x, doublereal *_y, doublereal alpha,
                 int *_le, doublereal *_q);

     void transform(doublereal *_xl, doublereal *_yl, doublereal *_zl,
                    doublereal *_xg, doublereal *_yg, doublereal *_zg,
                    doublereal *_str);

     void vonmis(doublereal rmx, doublereal rmy, doublereal rmxy,
                 doublereal rnx, doublereal rny, doublereal rnxy,
                 doublereal &t,   doublereal &sv);

     void compj2(doublereal sx, doublereal sy, doublereal sxy, doublereal &svm);
*/
};
 
#endif
