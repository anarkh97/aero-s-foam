#ifndef  _UTILITIES_H
#define  _UTILITIES_H

template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

void   orthonorm3( double rten[3][3] );
double fitalg3_3nodnew(double x0[3][3], double xn[3][3]);
void   inc_rottensor ( double[3], double[3][3]);
void   inc_rottensor ( double[3][3], double[3]);
double form_rottensor ( double[3], double[3][3]);
void   vec_to_mat (double rvec[3], double rten[3][3] );
void   vec_to_quat (double[3], double[4]);
void   quat_to_mat ( double[4], double[3][3]);
void   mat_to_vec( double[3][3], double[3]);
void   mat_to_quat( double[3][3], double[4]);
void   trans_fsl( double force[], double *stiff[],
               double t0n[3][3], int num_nodes );
void   pseudorot_var( double rvec[3], double varmat[3][3] );
void   pseudorot_2var( double r[3], double f[3], double scndvar[3][3]);
void   leftmult_rotvar(int  num_nodes, int itrans, double  rotvar[][3][3]
                                              , FullSquareMatrix &stiff );
void   rightmult_rotvar(int  num_nodes, int itran, double  rotvar[][3][3],
                                            FullSquareMatrix &stiff );
void   crossprod(double a[3], double b[3], double c[3]);
void   dcrossprod(double  a[3], double  b[3], double  c[3],
                  double da[3], double db[3], double dc[3]);
void   normalize(double a[3]);
void dnormalize(double v[3],double dv[3]);
double   magnitude(double a[3]);
void   tran_fsl(  double force[], FullSquareMatrix &stiff,
                double t0n[3][3], int num_nodes );
void   mat_mult_mat(const double R1[3][3], const double R2[3][3], double result[3][3], 
                    int transflag);
void mat_mult_vec(double A[3][3], double b[3], double c[3], int transflag = 0);

void tran_force(double* force, double tmat[3][3], int num_nodes,
                int num_dofs_per_node = 6 );
void tran_stiff(FullSquareMatrix &stiff, double tmat[3][3]);


#endif
