#ifndef _SHELLMATERIAL_HPP_
#define _SHELLMATERIAL_HPP_

#ifdef USE_EIGEN3
#include <iostream>
#include <vector>
#include <Eigen/Core>

template<typename doublereal>
class ShellMaterial
{
  public:
    virtual ~ShellMaterial() {}
    virtual void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                         doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0) = 0; 
                                                                          // Upsilon is the generalized "strains" {e,chi}
                                                                          // Sigma is the generalized "stress" {N,M}
                                                                          // D is the tangent constitutive matrix { Dm, Dmb; Dbm, Db }
                                                                          // temp is the temperature
    virtual void GetConstitutiveResponseSensitivityWRTthickness(doublereal *Upsilon, doublereal *Sigma, doublereal *dDdthick,
                                                                doublereal *eframe, int gp, doublereal temp = 0) { 
      std::cerr << "GetConstitutiveResponseSensitivityWRTthickness is not defined\n"; }
    virtual void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                           doublereal *eframe, int gp) {
      std::cerr << "GetConstitutiveResponseSensitivityWRTdisp is not defined\n"; }
    virtual doublereal* GetCoefOfConstitutiveLaw() { return NULL; }
    virtual doublereal GetShellThickness() = 0;
    virtual doublereal GetAreaDensity() = 0; // mass per unit area
    virtual doublereal GetAmbientTemperature() = 0;
    virtual void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                              doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0) = 0;
    virtual void GetLocalConstitutiveResponseSensitivityWRTthick(doublereal *_Upsilon, doublereal *_dsigmadh, doublereal dzdh,
                                                                 doublereal *, int) {
      std::cerr << "GetLocalConstitutiveResponseSensitivityWRTthick is not defined\n"; }
    virtual void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                                doublereal *eframe, int gp) {
      std::cerr << "GetLocalConstitutiveResponseSensitivityWRTdisp is not defined\n"; }
    virtual int GetNumStates() { return 0; }
    virtual void SetState(doublereal *state) {}
    virtual void GetState(doublereal *state) {}
    virtual void UpdateState(doublereal *Upsilon, doublereal *state, int gp, doublereal dt = 0) {}
    virtual std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z) { return std::vector<doublereal>(); }
    virtual std::vector<doublereal> GetLocalBackStress(int nd, doublereal z) { return std::vector<doublereal>(); }
    virtual doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z) { return 0; }
    virtual doublereal GetLocalDamage(int nd, doublereal z) { return 0; }
    virtual doublereal GetDissipatedEnergy(int gp) { return 0; }
    virtual bool CheckFailure() { return false; } // used to initiate element deletion

  protected:
    Eigen::Matrix<doublereal,3,3>
    andesinvt(doublereal *_eframe, doublereal *_aframe, doublereal thetaf);
};

//     ------------------------------------------------ 
//     ISOTROPIC (ESSENTIALLY PLANE) MATERIAL           
//     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS 
//     ------------------------------------------------ 
template<typename doublereal>
class ShellMaterialType0 : public ShellMaterial<doublereal>
{
    doublereal E, thick, nu;
    doublereal rho; // volume density
    doublereal Ta;  // ambient temperature
    doublereal w;   // coefficient of thermal expansion
  public:
    ShellMaterialType0(doublereal _E, doublereal _thick, doublereal _nu, doublereal _rho,
                       doublereal _Ta = 0., doublereal _w = 0.) 
      : E(_E), thick(_thick), nu(_nu), rho(_rho), Ta(_Ta), w(_w) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    void GetConstitutiveResponseSensitivityWRTthickness(doublereal *Upsilon, doublereal *Sigma, doublereal *dDdthick,
                                                        doublereal *eframe, int gp, doublereal temp = 0);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { 
      return rho*thick; 
    }
    doublereal GetAmbientTemperature() { return Ta; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    void GetLocalConstitutiveResponseSensitivityWRTthick(doublereal *_Upsilon, doublereal *_dsigmadh, doublereal dzdh,
                                                         doublereal *, int);
    void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                        doublereal *eframe, int gp);
};

//     ------------------------------------------------- 
//     COMPOSITE MATERIAL WITH KNOWN CONSTITUTIVE MATRIX 
//     ------------------------------------------------- 
template<typename doublereal>
class ShellMaterialType1 : public ShellMaterial<doublereal>
{
    Eigen::Map<Eigen::Matrix<doublereal,6,6,Eigen::RowMajor> > coef;
    doublereal *aframe;
    doublereal rhoh;
    doublereal thick;
    doublereal Ta; // ambient temperature
    Eigen::Map<Eigen::Matrix<doublereal,6,1> > Alpha; // coefficients of thermal expansion
    static bool Wlocal_stress;
  public:
    ShellMaterialType1(doublereal *_coef, doublereal *_aframe, doublereal _rhoh, doublereal _thick = 0.,
                       doublereal _Ta = 0.)
      : coef(_coef), aframe(_aframe), rhoh(_rhoh), thick(_thick), Ta(_Ta), Alpha(_coef+36) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    doublereal* GetCoefOfConstitutiveLaw() { return coef.data(); }
    doublereal GetShellThickness();
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetAmbientTemperature() { return Ta; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                              doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
};

//     ---------------------------------------------- 
//     COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES 
//     ---------------------------------------------- 
template<typename doublereal>
class ShellMaterialTypes2And3 : public ShellMaterial<doublereal>
{
    int nlayer;
    Eigen::Map<Eigen::Matrix<doublereal,12,Eigen::Dynamic> > mtlayer;
    bool couple;
    doublereal *aframe;
    doublereal thick;
    doublereal rhoh;
    doublereal Ta;
    doublereal nsm;

  public:
    ShellMaterialTypes2And3(int _nlayer, doublereal *_mtlayer, bool _couple, doublereal *_aframe, doublereal _Ta = 0.,
                            doublereal _nsm = 0.);

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetAmbientTemperature() { return Ta; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
};

//     ------------------------------------------------ 
//     ISOTROPIC LINEAR ELASTIC - J2 PLASTIC MATERIAL   
//     ------------------------------------------------ 
template<typename doublereal, typename localmaterial>
class ShellMaterialType4 : public ShellMaterial<doublereal>
{
    doublereal thick, nu;
    doublereal rho;
    localmaterial **mat;
    int nlayer; // number of material points through the thickness of the shell
    int maxgus; // number of material points over the area of the shell, per layer
    doublereal Ta; // ambient temperature
    doublereal w;  // coefficient of thermal expansion
  public:
    ShellMaterialType4(doublereal _thick, doublereal _nu, doublereal _rho, localmaterial *_mat, int _nlayer, int _maxgus,
                       doublereal _Ta = 0., doublereal _w = 0.)
      : thick(_thick), nu(_nu), rho(_rho), nlayer(_nlayer), maxgus(_maxgus), Ta(_Ta), w(_w) {
      mat = new localmaterial * [_nlayer*_maxgus];
      for (int i = 0; i < _nlayer*_maxgus; ++i) mat[i] = _mat->Clone();
    }
    ~ShellMaterialType4() {
      for(int i = 0; i < nlayer*maxgus; ++i) delete mat[i];
      delete [] mat;
    }

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { return rho*thick; }
    doublereal GetAmbientTemperature() { return Ta; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp, doublereal temp = 0, doublereal dt = 0);
    int GetNumStates() { return nlayer*maxgus*7; } // TODO 7 should be provided by the localmaterial
    void SetState(doublereal *state);
    void GetState(doublereal *state);
    void UpdateState(doublereal *Upsilon, doublereal *state, int gp, doublereal dt = 0);
    std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z);
    std::vector<doublereal> GetLocalBackStress(int nd, doublereal z);
    doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z);
    doublereal GetLocalDamage(int nd, doublereal z);
    doublereal GetDissipatedEnergy(int gp);
    bool CheckFailure();
};
#endif

#endif
