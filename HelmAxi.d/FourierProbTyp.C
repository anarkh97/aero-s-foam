#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <Driver.d/Domain.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Solvers.d/Solver.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/FourierDescrip.h>
#include <HelmAxi.d/FourierProbTyp.h>


FourierSolver::FourierSolver() {
}


FourierSolver::FourierSolver(FourierStatic *PrbD) {
   probDesc = PrbD;
   globalsol= NULL;
}


void FourierSolver::solve() {

 probDesc->preProcess();

 int i=0;
 int mode = probDesc->getModes();
 
 globalsol = new ComplexVector[2*mode+1];

 probDesc->constructK();

 fprintf(stderr," ... Running Fourier mode #0        ...\n");

 probDesc->BuildKs(0);

 ComplexSolver *mat = probDesc->getSolver();
 ComplexVector rhs(probDesc->solVecInfo(),DComplex(0,0));
 ComplexVector sol_n(probDesc->solVecInfo(),DComplex(0,0));

 mat->factor();

 probDesc->buildRHS(rhs,i);

 mat->solve(rhs,sol_n);

 globalsol[0]=probDesc->MergeSolBCs(sol_n);

 for (i = 1; i <= mode ; ++i) {

   fprintf(stderr," ... Running Fourier mode #%d        ...\n",i);

   probDesc->BuildKs(i);
   mat = probDesc->getSolver();
   mat->factor(); 

   rhs.zero();
   sol_n.zero();
   probDesc->buildRHS(rhs,2*i-1);

   mat->solve(rhs,sol_n);
   globalsol[2*i-1]=probDesc->MergeSolBCs(sol_n);

   rhs.zero();
   sol_n.zero();
   probDesc->buildRHS(rhs,2*i);
   mat->solve(rhs,sol_n);

   globalsol[2*i]=probDesc->MergeSolBCs(sol_n);

 }

 probDesc->postProcessing();
 probDesc->Reconstruction3D(globalsol);

/*
 // Output the Fourier coefficients of the solution
 //
 // UH - 12-17-99
 // Put this output file in a more general format

 int k,j;
 int maxNode = probDesc->getDomain()->numNode();

 FILE *fid = fopen("FourierCoef","w");

 for (k=0; k<(2*mode+1); ++k)
   for (j=0; j<maxNode; ++j) {
     fprintf(fid,"%1.7e",real(globalsol[k][j]));
     if (imag(globalsol[k][j]) < 0.0)
        fprintf(fid,"  %1.7e \n",imag(globalsol[k][j]));
     else
        fprintf(fid,"  %1.7e \n",imag(globalsol[k][j]));
   }

*/

}


