#ifndef _SHELLELEMENTTEMPLATE_HPP_
#define _SHELLELEMENTTEMPLATE_HPP_

#include <Element.d/FelippaShell.d/ShellMaterial.hpp>

// Note: the following needs to be defined for full compatability with the implementation of elements 8 and 20
//#define COMPATABILITY_MODE

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

    void
    andescrd(int elm, doublereal *x, doublereal *y, doublereal *z, doublereal *rot,
             doublereal *xlp, doublereal *ylp, doublereal *zlp, doublereal &area);

    void
    andesms(int elm, doublereal *x, doublereal *y, doublereal *z,
            doublereal *emass, doublereal *gamma, doublereal *grvfor,
            bool grvflg, doublereal &totmas, bool masflg);

    void
    andesstf(int elm, doublereal *estiff, doublereal *fint, doublereal nu,
             doublereal *x, doublereal *y, doublereal *z, doublereal *globalu,
             int ctyp, int flag);

    void
    andesvms(int elm, int maxstr, doublereal nu, doublereal *globalx,
             doublereal *globaly, doublereal *globalz, doublereal *globalu,
             doublereal *stress, int ctyp, int strainflg, int surface);

    doublereal
    equivstr(doublereal sxx, doublereal syy, doublereal szz, doublereal sxy);

    void
    transform(doublereal *lframe, doublereal *gframe, doublereal *str);

    void
    andesups(int elm, doublereal *state, doublereal *X, doublereal *Y, 
             doublereal *Z, doublereal *v);
};

#endif