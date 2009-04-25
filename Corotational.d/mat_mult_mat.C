
void
mat_mult_mat(double A[3][3], double B[3][3], double C[3][3], int transflag)
{
 int i,j;
 if(transflag == 0) {
   // C = A*B
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
 }
 else if(transflag == 1) {
   // C = A^t*B
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[0][i]*B[0][j] + A[1][i]*B[1][j] + A[2][i]*B[2][j];

 } 
 else if(transflag == 2) {
   // C = A*B^t
   for(i=0; i<3; ++i)
     for(j=0; j<3; ++j)
       C[i][j] = A[i][0]*B[j][0] + A[i][1]*B[j][1] + A[i][2]*B[j][2];

 }

}
