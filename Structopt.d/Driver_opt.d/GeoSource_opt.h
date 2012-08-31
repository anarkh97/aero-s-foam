#ifndef _GEO_SOURCE_OPT_HPP_
#define _GEO_SOURCE_OPT_HPP_

#ifdef STRUCTOPT

#include <map>
#include <Driver.d/GeoSource.h>


// Density projection stuff
struct DensityProjData
{
  int densProjFlag;
  double densProjParam1, densProjParam2, densProjParam3;
  DensityProjData() : densProjFlag(0), densProjParam1(0), densProjParam2(0), densProjParam3(0) {}
};

// Nodal density stuff
struct NodalDensityData
{
  int nodDensFunc;
  double nodDensParam1, nodDensParam2;
  NodalDensityData() : nodDensFunc(0), nodDensParam1(0), nodDensParam2(0) {}
};


class GeoSource_opt: public GeoSource 
{
private:
  SPropContainer sgradprops;   // set of gradients of properties
  bool movingNodes; // nodes change their positions??
  bool dsgnDepRHS;  // design-dep. RHS???
  std::map<int, CoefData> d_coefData;
  DensityProjData densProjData;
  NodalDensityData nodalDensData;

public:
  explicit GeoSource_opt(int iniSize = 16): 
    GeoSource(iniSize), movingNodes(true), dsgnDepRHS(true) {}
  virtual ~GeoSource_opt() {}

  SPropContainer& getGradStructProps() { return sgradprops; }
  CoefData* getGradCoefData(int i);
  void outputEnergies(int, double, double, double, double, double, double);

  // function for gradient arrays of structural properties
  void buildOptGradProp(); 
  void zeroGradProp();

  bool areThereMovingNodes() const { return movingNodes; }
  void setMovingNodesFlag(bool f) { movingNodes = f; return; }
  bool isRHSDsgnDep() const { return dsgnDepRHS; }
  void setDsgnDepRHSFlag(bool f) { dsgnDepRHS = f; return; }

  NodalDensityData& getNodalDensityData() { return nodalDensData; }
  DensityProjData& getDensityProjData() { return densProjData; }
};

#endif

#endif
