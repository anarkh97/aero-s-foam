// HB (06/11/05)
// To be included for solid element 
// EXPERIMENTAL ...
#include <Math.d/FullSquareMatrix.h>

// For stiffness & mass matrices 3D solid elements 
extern void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
extern void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);

extern void addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls);
extern void addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls);

extern int checkJacobian(double *J, int *jSign, char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);

// For ansitropic constitutive matrix
extern void rotateConstitutiveMatrix(double *_Cin, double *T33, double Cout[6][6]);

// For stresses & strains evaluation in case of ansitropic constitutive matrix
extern void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
extern double computeStress3DSolid(double Stress[6],double Strain[6], double C[6][6]);
extern double computeVonMisesStress(double Stress[6]);
extern double computeVonMisesStrain(double Strain[6]);
 
