// If the BlastLoading class is not defined, define it:
#ifndef _BLASTLOADING_H_
#define _BLASTLOADING_H_
class BlastLoading {
    public:
        struct BlastData {
            double ExplosivePosition[3];
            double ExplosiveDetonationTime;
            enum {SurfaceBurst, AirBurst} BlastType;
            double ExplosiveWeight;
            double ExplosiveWeightCubeRoot;
            double ScaleLength;
            double ScaleTime;
            double ScaleMass;
            bool ConwepGlobalOnOff;
            void print();
        };
    private:
        class Conwep {
            public:
                static double Blast(
                    const BlastLoading::BlastData& P,
                    const double CurrentElementFaceCentroidPosition[3],
                    const double CurrentElementFaceNormalDirection[3],
                    double CurrentTime);
                static double Decay(
                    double CurrentPressure,
                    double CurrentImpulse,
                    double PositivePhaseDuration);
                static double IncidentPressure(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ReflectedPressure(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ArrivalTime(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10) ;
                static double PositivePhaseDuration(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ReflectedImpulse(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double IncidentImpulse(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10) ;
                static void Parameters(
                    const BlastLoading::BlastData& P,
                    double DistanceFromElementFaceCentroidToExplosive,
                    double& IncidentWaveArrivalTime,
                    double& PositivePhaseDuration,
                    double& IncidentWaveImpulse,
                    double& ReflectedWaveImpulse,
                    double& IncidentWavePressure,
                    double& ReflectedWavePressure,
                    double& IncidentWaveDecayExponent,
                    double& ReflectedWaveDecayExponent);
                static double Pressure(
                    double CurrentTimeSinceExplosionTime,
                    double IncidentWaveArrivalTime,
                    double PositivePhaseDuration,
                    double IncidentWavePressure,
                    double ReflectedWavePressure,
                    double CurrentElementPositionCosine,
                    double IncidentWaveDecayExponent,
                    double ReflectedWaveDecayExponent);
        };
    public:
        static double ComputeShellPressureLoad(
            const double* CurrentElementNodePositions,
            double CurrentTime,
            const BlastLoading::BlastData& P );
        static BlastData InputFileData;
        static bool WarnedZeroDist;
        static bool WarnedDecayExp;
};
#endif
