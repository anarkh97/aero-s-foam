#ifndef COEFFOURIER_H_
#define COEFFOURIER_H_

#include<Utils.d/MyComplex.h>

DComplex coefExpDir(int mode, double k, double rho, double phi);

DComplex coefExpNeu(int mode, double k, double rho, double phi, double nr,
                    double nz, double dx, double dy, double dz);

double besselj(int n, double x);

#endif
