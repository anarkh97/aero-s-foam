#ifndef _FABRICMAT_H_
#define _FABRICMAT_H_

#include <Element.d/NonLinearity.d/NLMaterial.h>

class StructProp;

class FabricMat : public NLMaterial
{
protected:
	int pxx_map_id, pyy_map_id, sxy_map_id;
	double rho, t, Tref, alpha;
	SS2DTData *pxx_map, *pyy_map;
	MFTTData *sxy_map;

public:
	FabricMat(StructProp *p);
	FabricMat(double _rho, int _pxx_map_id, int _pyy_map_id, int _sxy_map_id, double _t, double _Tref, double _alpha);

	int getNumStates() const override { return 0; }

	void getStress(Tensor *stress, Tensor &strain, double*, double temp) override;

	void getTangentMaterial(Tensor *tm, Tensor &strain, double*, double temp) override;

	void getElasticity(Tensor *tm) const override {};

	void updateStates(Tensor &en, Tensor &enp, double *state, double temp) override {};

	void getStressAndTangentMaterial(Tensor *stress, Tensor *tm, Tensor &strain, double*, double temp) override;

	void integrate(Tensor *stress, Tensor *tm, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp,
	               Tensor *cache, double dt=0) const override;

	void integrate(Tensor *stress, Tensor &en, Tensor &enp,
	               double *staten, double *statenp, double temp,
	               Tensor *cache, double dt=0) const override;

	void initStates(double *) override {};

	GenStrainEvaluator<TwoDTensorTypes<9> > * getGenStrainEvaluator();

	double getDensity() override { return rho; }

	double getThickness() { return t; }

	double getReferenceTemperature() override { return Tref; }

	void setS1DProps(MFTTData *ss1dt) override;
	void setS2DProps(SS2DTData *ss2dt) override;
};

#endif
