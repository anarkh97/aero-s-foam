// BlastLoading.C
// This file implements the CONWEP blast loading model.
#include "BlastLoading.h"
#include <math.h>
#include <iostream>
#include <fstream>
// Variables and Notes:
// P: contains the input file parameters: explosive location, detonation time, explosion type (air blast or ground blast), charge weight, and the cube root of the charge weight.
// R: distance of current point to the explosive location.
// arrivalTime: time for the shock wave to arrive at the surface.
// positivePhaseDuration: time during which the shock wave is above ambient pressure.
// There are 2 shock waves: the incident normal shock wave, and the ground-reflected oblique shock wave.
// incidentImpulse and incidentPressure: specific impulse and pressure of the incident shock wave.
// reflectedImpulse and reflectedPressure: specific impulse and pressure of the ground-reflected shock wave.
// a: decay exponent of the incident shock wave.
// b: decay exponent of the ground-reflected shock wave.


// Define the BlastLoading::Conwep::Blast function, which returns a pressure:
// P contains everything in BlastData: x0, t0, blastType, chargeWeight, chargeWeightCubeRoot, scaleLength, scaleTime, scaleMass.
double BlastLoading::Conwep::Blast(const BlastLoading::BlastData& P,
                                   const double x[3], // Element face centroid
                                   const double n[3], // Element face normal
                                   double t) { // Current time.
  double f[3] = {
    P.x0[0]-x[0],
    P.x0[1]-x[1],
    P.x0[2]-x[2]
  };
  double R = sqrt( f[0]*f[0]+f[1]*f[1]+f[2]*f[2] );
  f[0] /= R;
  f[1] /= R;
  f[2] /= R;
  // Convert distance to feet:
  //R *= 3.2808399;
  R = R/P.scaleLength*3.2808399;
  double posCosine = n[0]*f[0]+n[1]*f[1]+n[2]*f[2];
  double arrivalTime;
  double positivePhaseDuration;
  double incidentImpulse;
  double reflectedImpulse;
  double incidentPressure;
  double reflectedPressure;
  double a;
  double b;
  Conwep::Params(P,
                 R,
                 arrivalTime,
                 positivePhaseDuration,
                 incidentImpulse,
                 reflectedImpulse,
                 incidentPressure,
                 reflectedPressure,
                 a,
                 b);
  //double ts = (t- P.t0)*1000.0; // Convert time to milliseconds.
  double ts = (t- P.t0)*1000.0*P.scaleTime;
  double p = Conwep::Pressure(ts,
                              arrivalTime,
                              positivePhaseDuration,
                              incidentPressure,
                              reflectedPressure,
                              posCosine,
                              a,
                              b);
  return p;
}
double BlastLoading::Conwep::Decay(double p0,
                                   double i0,
                                   double td) {
// p0 refers to either incidentPressure or reflectedPressure, depending on which one called the function.
// i0 refers to either incidentImpulse or reflectedImpulse, depending on which one called the function.
// td refers to positivePhaseDuration.
  double ptoi = p0*td/i0;
  double a = ptoi-1.0;
  double fa,fpa;
  do {
    fa = a*a-ptoi*(a+exp(-a)-1.0);
    fpa = 2.0*a-ptoi*(1.0-exp(-a));
    a = a-fa/fpa;
  } while (fabs(fa) > 1.0e-6);
  return a;
}
double BlastLoading::Conwep::IncidentPressure(const BlastLoading::BlastData& P,
                                              double zlog) {
  const static double cpso[2][12] = {  1.9422502013,
                                      -1.6958988741,
                                      -0.154159376846,
                                       0.514060730593,
                                       0.0988534365274,
                                      -0.293912623038,
                                      -0.0268112345019,
                                       0.109097496421,
                                       0.00162846756311,
                                      -0.0214631030242,
                                       0.0001456723382,
                                       0.00167847752266,
                                       1.77284970457,
                                      -1.69012801396,
                                       0.00804973591951,
                                       0.336743114941,
                                      -0.00516226351334,
                                      -0.0809228619888,
                                      -0.00478507266747,
                                       0.00793030472242,
                                       0.0007684469735,
                                       0.0,
                                       0.0,
                                       0.0 };
  double u = -0.756579301809 + 1.35034249993*zlog;
  int i = (int)P.blastType;
  double pinc = cpso[i][11];
  for (int j = 10; j >= 0; --j) {
    pinc = pinc*u+cpso[i][j];
  }
  return pow(10.0, pinc);
}
double BlastLoading::Conwep::ReflectedPressure(const BlastLoading::BlastData& P,
					       double zlog) {
  const static double csurf[12] = {  2.56431321138,
                                    -2.21030870597,
                                    -0.218536586295,
                                     0.895319589372,
                                     0.24989009775,
                                    -0.569249436807,
                                    -0.11791682383,
                                     0.224131161411,
                                     0.0245620259375,
                                    -0.0455116002694,
                                    -0.00190930738887,
                                     0.00361471193389 };
  const static double cfree[10] = {  2.39106134946,
                                    -2.21400538997,
                                     0.035119031446,
                                     0.657599992109,
                                     0.0141818951887,
                                    -0.243076636231,
                                    -0.0158699803158,
                                     0.0492741184234,
                                     0.00227639644004,
                                    -0.00397126276058 };
  double u;
  double pref;
  int j;
  switch (P.blastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    u = -0.789312405513+1.36637719229*zlog;
    pref = csurf[11];
    for (j = 10; j >= 0; --j)
      pref = pref*u+csurf[j];
    break;
  case BlastLoading::BlastData::AirBurst:
    u = -0.756579301809 + 1.35034249993*zlog;
    pref = cfree[9];
    for (j = 8; j >= 0; --j)
      pref = pref*u+cfree[j];
    break;  
  }
  return pow(10.0, pref);
}
double BlastLoading::Conwep::ArrivalTime(const BlastLoading::BlastData& P,
					 double zlog) {
  const static double csurf[10] = { -0.173607601251,
                                     1.35706496258,
                                     0.052492798645,
                                    -0.196563954086,
                                    -0.0601770052288,
                                     0.0696360270891,
                                     0.0215297490092,
                                    -0.0161658930785,
                                    -0.00232531970294,
                                     0.00147752067524 };
  const static double cfree[8] = { -0.0423733936826,
                                    1.36456871214, 
                                   -0.0570035692784,
                                   -0.182832224796,
                                    0.0118851436014,
                                    0.0432648687627,
                                   -0.0007997367834,
                                   -0.00436073555033 };
  double u;
  double tarr;
  int j;
  switch (P.blastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    u = -0.755684472698 + 1.37784223635*zlog;
    tarr = csurf[9];
    for (j = 8; j >= 0; --j)
      tarr = tarr*u+csurf[j];
    break;
  case BlastLoading::BlastData::AirBurst:
    u = -0.80501734056 + 1.37407043777*zlog;
    tarr = cfree[7];
    for (j = 6; j >= 0; --j)
      tarr = tarr*u+cfree[j];
    break;  
  }
  return pow(10.0, tarr);
}
double BlastLoading::Conwep::PositivePhaseDuration(const BlastLoading::BlastData& P,
                                                   double zlog) {
  const static double csurf1[6] = { -0.728671776005,
                                     0.130143717675,
                                     0.134872511954,
                                     0.0391574276906,
                                    -0.00475933664702,
                                    -0.00428144598008 };
  const static double csurf2[9] = {  0.20096507334,
                                    -0.0297944268976,
                                     0.030632954288,
                                     0.0183405574086,
                                    -0.0173964666211,
                                    -0.00106321963633,
                                     0.00562060030977,
                                     0.0001618217499,
                                    -0.0006860188944 };
  const static double csurf3[6] = {  0.572462469964,
                                     0.0933035304009,
                                    -0.0005849420883,
                                    -0.00226884995013,
                                    -0.00295908591505,
                                     0.00148029868929 };
  const static double cfree1[9] = { -0.801052722864,
                                     0.164953518069,
                                     0.127788499497,
                                     0.00291430135946,
                                     0.00187957449227,
                                     0.0173413962543,
                                     0.00269739758043,
                                    -0.00361976502798,
                                    -0.00100926577934 };
  const static double cfree2[9] = {  0.115874238335,
                                    -0.0297944268969,
                                     0.0306329542941,
                                     0.018340557407,
                                    -0.0173964666286,
                                    -0.00106321963576,
                                     0.0056206003128,
                                     0.0001618217499,
                                    -0.0006860188944 };
  const static double cfree3[8] = {  0.50659210403,
                                     0.0967031995552,
                                    -0.00801302059667,
                                     0.00482705779732,
                                     0.00187587272287,
                                    -0.00246738509321,
                                    -0.000841116668,
                                     0.0006193291052 };
  double u;
  double tdur;
  int j;
  switch (P.blastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    if (zlog <= -0.34) {
      tdur = -0.725;
    } else if (zlog <= 0.4048337) {
      u = -0.1790217052 + 5.25099193925*zlog;
      tdur = csurf1[5];
      for (j = 4; j >= 0; --j)
        tdur = tdur*u+csurf1[j];
    } else if (zlog <= 0.845098) {
      u = -5.85909812338 + 9.2996288611*zlog;
      tdur = csurf2[8];
      for (j = 7; j >= 0; --j)
        tdur = tdur*u+csurf2[j];
    } else  {
      u = -4.92699491141 + 3.46349745571*zlog;
      tdur = csurf3[5];
      for (j = 4; j >= 0; --j)
        tdur = tdur*u+csurf3[j];
    }
    break;
  case BlastLoading::BlastData::AirBurst:  
    if (zlog <= -0.34) {
      tdur = -0.824;
    } else if (zlog <= 0.350248) {
      u = 0.209440059933 + 5.11588554305*zlog;
      tdur = cfree1[8];
      for (j = 7; j >= 0; --j)
        tdur = tdur*u+cfree1[j];
    } else if (zlog <= 0.7596678) {
      u = -5.06778493835 + 9.2996288611*zlog;
      tdur = cfree2[8];
      for (j = 7; j >= 0; --j)
        tdur = tdur*u+cfree2[j];
    } else  {
      u = -4.39590184126 + 3.1524725264*zlog;
      tdur = cfree3[7];
      for (j = 6; j >= 0; --j)
        tdur = tdur*u+cfree3[j];
    }
    break;
  }
  return pow(10.0, tdur);
}
double BlastLoading::Conwep::ReflectedImpulse(const BlastLoading::BlastData& P,
					      double zlog) {
  const static double csurf[4] = { 1.75291677799,
                                  -0.949516092853,
                                   0.112136118689,
                                  -0.0250659183287 };
  const static double cfree[4] = { 1.60579280091,
                                  -0.903118886091,
                                   0.101771877942,
                                  -0.0242139751146 };
  double u;
  double ximpr;
  int j;
  switch (P.blastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    u = -0.781951689212 + 1.33422049854*zlog;
    ximpr = csurf[3];
    for (j = 2; j >= 0; --j)
      ximpr = ximpr*u+csurf[j];
    break;
  case BlastLoading::BlastData::AirBurst:
    u = -0.757659920369 + 1.37882996018*zlog;
    ximpr = cfree[3];
    for (j = 2; j >= 0; --j)
      ximpr = ximpr*u+cfree[j];
    break;  
  }
  return pow(10.0, ximpr);
}
double BlastLoading::Conwep::IncidentImpulse(const BlastLoading::BlastData& P,
                                             double zlog) {
  const static double csurf1[5] = { 1.57159240621,
                                   -0.502992763686,
                                    0.171335645235,
                                    0.0450176963051,
				   -0.0118964626402 };
  const static double csurf2[8] = { 0.719852655584,
                                   -0.384519026965,
                                   -0.0280816706301,
                                    0.00595798753822,
                                    0.014544526107,
                                   -0.00663289334734,
                                   -0.00284189327204,
                                    0.0013644816227 };
  const static double cfree1[5] = { 1.43534136453,
                                   -0.443749377691,
                                    0.168825414684,
                                    0.0348138030308,
                                   -0.010435192824 };
  const static double cfree2[9] = { 0.599008468099,
                                   -0.40463292088,
                                   -0.0142721946082,
                                    0.00912366316617,
                                   -0.0006750681404,
                                   -0.00800863718901,
                                    0.00314819515931,
                                    0.00152044783382,
                                   -0.0007470265899};
  double u;
  double ximps;
  int j;
  switch (P.blastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    if (zlog <= 0.382017) {
      u = 0.832468843425 + 3.0760329666*zlog;
      ximps = csurf1[4];
      for (j = 3; j >= 0; --j)
        ximps = ximps*u+csurf1[j];
    } else  {
      u = -2.91358616806 + 2.40697745406*zlog;
      ximps = csurf2[7];
      for (j = 6; j >= 0; --j)
        ximps = ximps*u+csurf2[j];
    }
    break;
  case BlastLoading::BlastData::AirBurst:  
    if (zlog <= 0.30103) {
      u = 1.04504877747 + 3.24299066475*zlog;
      ximps = cfree1[4];
      for (j = 3; j >= 0; --j)
        ximps = ximps*u+cfree1[j];
    } else  {
      u = -2.67912519532 + 2.30629231803*zlog;
      ximps = cfree2[8];
      for (j = 7; j >= 0; --j)
        ximps = ximps*u+cfree2[j];
    }
    break;
  }
  return pow(10.0, ximps);
}
void BlastLoading::Conwep::Params(const BlastLoading::BlastData& P,
                                  double R,
                                  double& arrivalTime,
                                  double& positivePhaseDuration,
                                  double& incidentImpulse,
                                  double& reflectedImpulse,
                                  double& incidentPressure,
                                  double& reflectedPressure,
                                  double& a,
                                  double& b) {
  static int cnt = 0;
  double z = R / P.chargeWeightCubeRoot;
  static std::ofstream det("DetonationProperties.txt");
  double zlog = log10(z);
  double zlo = (P.blastType == BlastLoading::BlastData::SurfaceBurst ? 0.45 : 0.37);
  arrivalTime = Conwep::ArrivalTime(P,zlog) * P.chargeWeightCubeRoot;
  positivePhaseDuration = Conwep::PositivePhaseDuration(P,zlog) * P.chargeWeightCubeRoot;
  incidentImpulse = Conwep::IncidentImpulse(P, zlog) * P.chargeWeightCubeRoot;
  reflectedImpulse = Conwep::ReflectedImpulse(P, zlog) * P.chargeWeightCubeRoot;
  incidentPressure = Conwep::IncidentPressure(P, zlog);
  reflectedPressure = Conwep::ReflectedPressure(P, zlog);
  //if ((++cnt) < 122) {
    det << "R = " << R << "\n"
        << "Arrival Time = " << arrivalTime << "\n"
        << "Positive Phase Duration = " << positivePhaseDuration << "\n"
        << "Incident Impulse = " << incidentImpulse << "\n"
        << "Reflected Impulse = " << reflectedImpulse << "\n"
        << "Incident Pressure = " << incidentPressure << "\n"
        << "Reflected Pressure = " << reflectedPressure << std::endl;
  //}
  if (z >= zlo) {
    a = Conwep::Decay(incidentPressure, incidentImpulse, positivePhaseDuration);
    b = Conwep::Decay(reflectedPressure, reflectedImpulse, positivePhaseDuration); 
  } else {
    a = b = 0;
  }
}
double BlastLoading::Conwep::Pressure(double ts,
                                      double arrivalTime,
                                      double positivePhaseDuration,
                                      double incidentPressure,
                                      double reflectedPressure,
                                      double posCosine,
                                      double a,
                                      double b) {
  if (ts >= arrivalTime) {
    double exa = exp(-a*(ts-arrivalTime)/positivePhaseDuration);
    double exb = exp(-b*(ts-arrivalTime)/positivePhaseDuration);
    double poscosa = (posCosine>0.0?posCosine:0.0);
    double p = (incidentPressure*exa*(1.0+poscosa-2.0*poscosa*poscosa)+reflectedPressure*exb*poscosa*poscosa)*(1.0-(ts-arrivalTime)/positivePhaseDuration);
    return (p>-14.7?p:-14.7);
  } else
    return 0.0;  
}
double BlastLoading::ComputeShellPressureLoad(const double* coords,
                                              double currentTime,
                                              const BlastLoading::BlastData& P ) {
  double a[3] = {
    coords[6]-coords[0],
    coords[7]-coords[1],
    coords[8]-coords[2]
  };
  double b[3] = {
    coords[9]-coords[3],
    coords[10]-coords[4],
    coords[11]-coords[5]
  };
  double n[3] = {
    a[1]*b[2]-a[2]*b[1],
    a[2]*b[0]-a[0]*b[2],
    a[0]*b[1]-a[1]*b[0]
  };
  double magn = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] /= magn;
  n[1] /= magn;
  n[2] /= magn;
  double x[3] = {0,0,0};
  for (int k = 0; k < 4; ++k) { 
    for (int j = 0; j < 3; ++j)
      x[j] += coords[k*3+j];
  }
  for (int j = 0; j < 3; ++j)
    x[j] *= 0.25;
  double p = Conwep::Blast(myData,x,n,currentTime);
  // p is in psi: convert it to Pa, then use scaleLength, scaleTime and scaleMass to convert it to correct pressure units:
  //return -p*6.895e3 // Convert psi to Pa.
  return -p*6.8947573e3/P.scaleMass*P.scaleLength*P.scaleTime*P.scaleTime;
}
// Initialize myData:
BlastLoading::BlastData BlastLoading::myData = {{0.0,0.0,0.0},0.0,
                                                BlastLoading::BlastData::AirBurst,1.0,0.0,0.3048,1.0,1.0};
