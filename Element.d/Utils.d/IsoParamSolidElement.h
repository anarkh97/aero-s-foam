// --------------------------------------------------------------------
// HB - 05/26/05
// --------------------------------------------------------------------
// Temptative of a base class for (mainly 3D) isoparametric solid 
// elements.
// This can be usefull because some methods like computing the thermal 
// forces or the nodal stresses or strains can be nearly "blind" of the
// underlaying type of the element: just the integration rule & the
// shape fcts (& derivative) actually depends on the underlaying type 
// of the element.
// --------------------------------------------------------------------
// EXPERIMENTAL ...
// --------------------------------------------------------------------

#ifndef _ISOPARAM_SOLID_ELEMENT_H_
#define _ISOPARAM_SOLID_ELEMENT_H_

#include <Element.d/Element.h>

class IsoParamSolidElement:
{
  public:
    virtual void getShapeFct(double* Shape, double* m);
    virtual void getdShapeFct(double (*dShape)[3], double* m);
    virtual double getDShapeFct(double (*dShape)[3], double* m, CoordSet& cs, double (*DShape)[3]=0);
    virtual double getDShapeFct(double (*dShape)[3], double* m, double* X, double* Y, double* Z, double (*DShape)[3]=0);
    virtual double getJacobian(double (*dShape)[3], double* m, CoordSet &cs, double (*DShape)[3]=0);
     
    virtual int getNodalStress();
    virtual int getNodalStrain();
    
    virtual int getThermalForces(Vector& elThermalForces,Vector& ndTemps, CoordSet& cs);
}
#endif

