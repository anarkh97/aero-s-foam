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
                                         doublereal *eframe, int gp) = 0; // Upsilon is the generalized "strains" {e,chi}
                                                                          // Sigma is the generalized "stress" {N,M}
                                                                          // D is the tangent constitutive matrix { Dm, Dmb; Dbm, Db }
    virtual void GetConstitutiveResponseSensitivityWRTthickness(doublereal *Upsilon, doublereal *Sigma, doublereal *dDdthick,
                                                                doublereal *eframe, int gp) { 
      std::cerr << "GetConstitutiveResponseSensitivityWRTthickness is not defined\n"; }
    virtual void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                           doublereal *eframe, int gp) {
      std::cerr << "GetConstitutiveResponseSensitivityWRTdisp is not defined\n"; }
    virtual void GetLocalConstitutiveResponseSensitivityWRTthick(doublereal *_Upsilon, doublereal *_dsigmadh, doublereal dzdh,
                                                                 doublereal *, int) {
      std::cerr << "GetLocalConstitutiveResponseSensitivityWRTthick is not defined\n"; }
    virtual void setThickness(doublereal _thick) = 0;
    virtual doublereal* GetCoefOfConstitutiveLaw() { 
      std::cerr << "GetCoefOfConstitutiveLaw() is not defined\n";  
      Eigen::Matrix<doublereal,36,1> A;
      A.setZero();
      return A.data(); 
    }
    virtual void resetLayerThickness(doublereal *layerthickness) { }
    virtual void resetAreaDensity() { }
    virtual void GetLayerThickness(doublereal *layerThickness) { }
    virtual doublereal GetShellThickness() = 0;
    virtual doublereal GetAreaDensity() = 0; // mass per unit area
    virtual doublereal GetSumDensity() { return 0; }
    virtual void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                              doublereal *eframe, int gp)
      { std::cerr << "GetLocalConstitutiveResponse is not defined\n"; }
    virtual void GetLocalConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dsigmadu, doublereal z,
                                                                doublereal *eframe, int gp)
      { std::cerr << "GetLocalConstitutiveResponseWRTdisp is not defined\n"; }
    virtual int GetNumStates() { return 0; }
    virtual void SetState(doublereal *state) {}
    virtual void GetState(doublereal *state) {}
    virtual void UpdateState(doublereal *Upsilon, doublereal *state, int gp) {}
    virtual std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z) { return std::vector<doublereal>(); }
    virtual std::vector<doublereal> GetLocalBackStress(int nd, doublereal z) { return std::vector<doublereal>(); }
    virtual doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z) { return 0; }
    virtual doublereal GetDissipatedEnergy(int gp) { return 0; }

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
  public:
    ShellMaterialType0(doublereal _E, doublereal _thick, doublereal _nu, doublereal _rho) 
      : E(_E), thick(_thick), nu(_nu), rho(_rho) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp);
    void GetConstitutiveResponseSensitivityWRTthickness(doublereal *Upsilon, doublereal *Sigma, doublereal *dDdthick,
                                                        doublereal *eframe, int gp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void GetLocalConstitutiveResponseSensitivityWRTthick(doublereal *_Upsilon, doublereal *_dsigmadh, doublereal dzdh,
                                                         doublereal *, int);
    void setThickness(doublereal _thick) { thick = _thick; }
    doublereal* GetCoefOfConstitutiveLaw() { 
      Eigen::Matrix<doublereal,36,1> A;
      A.setZero();
      return A.data(); 
    }
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { 
//      std::cerr << "rho is " << rho << std::endl;
//      std::cerr << "thick is " << thick << std::endl;
      return rho*thick; 
    }
    doublereal GetSumDensity() { return rho; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp);
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
  public:
    ShellMaterialType1(doublereal *_coef, doublereal *_aframe, doublereal _rhoh, doublereal _thick = 0.)
      : coef(_coef), aframe(_aframe), rhoh(_rhoh), thick(_thick) {}

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void setThickness(doublereal _thick) { thick = _thick; }
    doublereal* GetCoefOfConstitutiveLaw() { return coef.data(); }
    doublereal GetShellThickness();
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetSumDensity() { return rhoh/thick; } 
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
    doublereal rho;

  public:
    ShellMaterialTypes2And3(int _nlayer, doublereal *_mtlayer, bool _couple, doublereal *_aframe);

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp);
    void GetConstitutiveResponseSensitivityWRTdisp(doublereal *dUpsilondu, doublereal *dSigmadu, doublereal *D,
                                                   doublereal *eframe, int gp);
    void setThickness(doublereal _thick) { thick = _thick; }
    doublereal* GetCoefOfConstitutiveLaw() { 
      Eigen::Matrix<doublereal,36,1> A;
      A.setZero();
      return A.data(); 
    }
    void resetLayerThickness(doublereal *layerthickness) { for(int i=0; i<nlayer; ++i) mtlayer(7,i) = layerthickness[i]; }
    void resetAreaDensity() { rhoh = 0; for(int i=0; i<nlayer; ++i) rhoh += mtlayer(6, i)*mtlayer(7, i); }
    void GetLayerThickness(doublereal *layerThickness) { for(int i=0; i<nlayer; ++i) layerThickness[i] = mtlayer(7,i); }
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { return rhoh; }
    doublereal GetSumDensity() { return rho; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp);
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
  public:
    ShellMaterialType4(doublereal _thick, doublereal _nu, doublereal _rho, localmaterial *_mat, int _nlayer, int _maxgus)
      : thick(_thick), nu(_nu), rho(_rho), nlayer(_nlayer), maxgus(_maxgus) {
      mat = new localmaterial * [_nlayer*_maxgus];
      for (int i = 0; i < _nlayer*_maxgus; ++i) mat[i] = _mat->Clone();
    }
    ~ShellMaterialType4() {
      for(int i = 0; i < nlayer*maxgus; ++i) delete mat[i];
      delete [] mat;
    }

    void GetConstitutiveResponse(doublereal *Upsilon, doublereal *Sigma, doublereal *D,
                                 doublereal *eframe, int gp);
    void setThickness(doublereal _thick) { thick = _thick; }
    doublereal* GetCoefOfConstitutiveLaw() { 
      Eigen::Matrix<doublereal,36,1> A;
      A.setZero();
      return A.data(); 
    }
    doublereal GetShellThickness() { return thick; }
    doublereal GetAreaDensity() { return rho*thick; }
    doublereal GetSumDensity() { return rho; }
    void GetLocalConstitutiveResponse(doublereal *Upsilon, doublereal *sigma, doublereal z,
                                      doublereal *eframe, int gp);
    int GetNumStates() { return nlayer*maxgus*7; } // TODO 7 should be provided by the localmaterial
    void SetState(doublereal *state);
    void GetState(doublereal *state);
    void UpdateState(doublereal *Upsilon, doublereal *state, int gp);
    std::vector<doublereal> GetLocalPlasticStrain(int nd, doublereal z);
    std::vector<doublereal> GetLocalBackStress(int nd, doublereal z);
    doublereal GetLocalEquivalentPlasticStrain(int nd, doublereal z);
    doublereal GetDissipatedEnergy(int gp);
};
#endif

#endif
