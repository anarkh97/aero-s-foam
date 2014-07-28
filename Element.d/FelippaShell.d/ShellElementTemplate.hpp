#ifndef _SHELLELEMENTTEMPLATE_HPP_
#define _SHELLELEMENTTEMPLATE_HPP_

#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

// Note: The following needs to be defined for full compatability with the implementation of elements 8 and 20.
//       The recommended setting is undefined. In this case the stress output is more accurate, although a little
//       more costly, and the augmented mass matrix is used by default resulting in a larger critical time-step
//       for explicit dynamics.
#define COMPATIBILITY_MODE

template<typename doublereal, template<typename> class Membrane, template<typename> class Bending>
class ShellElementTemplate : public Membrane<doublereal>, public Bending<doublereal>
{
  protected:
    ShellMaterial<doublereal> *gpmat; // gauss points material
    ShellMaterial<doublereal> *nmat;  // nodes material

  public:
    ShellElementTemplate() : gpmat(NULL), nmat(NULL) {}

    ~ShellElementTemplate() {
      if(nmat && nmat != gpmat) delete nmat;
      if(gpmat) delete gpmat;
    }

    void setnmat(ShellMaterial<doublereal> *_mat) { nmat = _mat; }
    void setgpmat(ShellMaterial<doublereal> *_mat) { gpmat = _mat; }
    void setgpnmat(ShellMaterial<doublereal> *_mat) { 
      if(gpmat) delete gpmat;
      nmat = gpmat = _mat; 
    }

    void
    andescrd(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *rot,
             doublereal *xlp, doublereal *ylp, doublereal *zlp, doublereal &area);

    void
    andescrdWRTcoord(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *rot,
                     doublereal *xlp, doublereal *ylp, doublereal *zlp, doublereal &area,
                     doublereal *_dxpdx, doublereal *_dypdx, doublereal *_dzpdx,
                     doublereal *_dxlpdx, doublereal *_dylpdx, doublereal *_dzlpdx, doublereal *_dareadx);

    void
    andesgf(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *_gravityForce,
            doublereal *emass, doublereal *gamma, doublereal *grvfor,
            bool grvflg, doublereal &totmas, bool masflg, int gravflg);

    void
    andesms(int elm, doublereal *x, doublereal *y, doublereal *z,
            doublereal *emass, doublereal *gamma, doublereal *grvfor,
            bool grvflg, doublereal &totmas, bool masflg);

    void
    andesmsWRTthic(int elm, doublereal *x, doublereal *y, doublereal *z,
                   doublereal *gamma, doublereal *grvforSen, bool grvflg,
                   doublereal &totmasSen, bool masflg);

    void
    andesstf(int elm, doublereal *estiff, doublereal *fint, doublereal nu,
             doublereal *x, doublereal *y, doublereal *z, doublereal *globalu,
             int ctyp, int flag, doublereal *ndtemps = 0);

    void 
    andesstfWRTthick(int elm, doublereal *_destiffdthick, doublereal nu,
                     doublereal *x, doublereal *y, doublereal *z,
                     int ctyp, int flag);


    void
    andesvms(int elm, int maxstr, doublereal nu, doublereal *globalx,
             doublereal *globaly, doublereal *globalz, doublereal *globalu,
             doublereal *stress, int ctyp, int strainflg, int surface,
             doublereal *ndtemps = 0);

    void
    andesvmsWRTcoord(int elm, int maxstr, doublereal nu, doublereal *globalx,
                     doublereal *globaly, doublereal *globalz, doublereal *globalu,
                     doublereal *stress, doublereal *_dstressdx, int ctyp, int strainflg, int surface);

    void
    andesvmsWRTdisp(int elm, int maxstr, doublereal nu, doublereal *globalx,
             doublereal *globaly, doublereal *globalz, doublereal *globalu,
             doublereal *stress, doublereal *_vmsWRTdisp, int ctyp, int strainflg, int surface);

    void
    andesvmsWRTthic(int elm, int maxstr, doublereal nu, doublereal *globalx,
                    doublereal *globaly, doublereal *globalz, doublereal *globalu,
                    doublereal *stress, doublereal *_vmsWRTthic, int ctyp, int strainflg, int surface);

    doublereal
    equivstr(doublereal sxx, doublereal syy, doublereal szz, doublereal sxy);

    Eigen::Matrix<doublereal,1,18>
    equivstrSensitivityWRTdisp(doublereal vms, doublereal sxx, doublereal syy, doublereal szz, doublereal sxy,
                               Eigen::Matrix<doublereal,3,18> dsigmadu);

    Eigen::Matrix<doublereal,1,1>
    equivstrSensitivityWRTthic(doublereal vms, doublereal sxx, doublereal syy, doublereal szz, doublereal sxy,
                               Eigen::Matrix<doublereal,3,1> dsigmadh);

    void
    transform(doublereal *lframe, doublereal *gframe, doublereal *str);

    void
    andesups(int elm, doublereal *state, doublereal *X, doublereal *Y, 
             doublereal *Z, doublereal *v);

    void
    andesare(int elm, doublereal *X, doublereal *Y, doublereal *Z, doublereal &area);

    void
    andesden(int elm, doublereal *X, doublereal *Y, doublereal *Z, doublereal &D);
};

#endif
