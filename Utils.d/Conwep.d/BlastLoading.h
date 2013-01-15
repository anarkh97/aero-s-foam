// If the BlastLoading class is not defined, define it:
#ifndef _BLASTLOADING_H_
#define _BLASTLOADING_H_
class BlastLoading {
  public:
  struct BlastData {
    double x0[3];
    double t0;
    enum {SurfaceBurst, AirBurst} blastType;
    double chargeWeight,chargeWeightCubeRoot;
  };
  private:
  class Conwep {
  public:
    static double Blast(const BlastLoading::BlastData& P,
			const double x[3], // face centroid
			const double n[3], // face normal
			double t);
    static double Decay(double p0, double i0, double td);
    static double IncidentPressure(const BlastLoading::BlastData& P,
				   double zlog);
    static double ReflectedPressure(const BlastLoading::BlastData& P,
				    double zlog);
    static double ArrivalTime(const BlastLoading::BlastData& P,
			      double zlog) ;
    static double PositivePhaseDuration(const BlastLoading::BlastData& P,
					double zlog);
    static double ReflectedImpulse(const BlastLoading::BlastData& P,
				   double zlog);
    static double IncidentImpulse(const BlastLoading::BlastData& P,
				  double zlog) ;
    static void Params(const BlastLoading::BlastData& P,
		       double R,
		       double& arrivalTime,
		       double& positivePhaseDuration,
		       double& incidentImpulse,
		       double& reflectedImpulse,
		       double& incidentPressure,
		       double& reflectedPressure,
		       double& a, double& b);
    static double Pressure(double ts,double arrivalTime,
			   double positivePhaseDuration,
			   double incidentPressure,
			   double reflectedPressure,
			   double posCosine,
			   double a,double b);
  };
  public:
  static double ComputeShellPressureLoad(const double* coords, double currentTime );
  static BlastData myData;
//  static double currentTime;
};
#endif // _BLASTLOADING_H_
