#include <Structopt.d/Driver_opt.d/StaticProbType_opt.h>

#ifdef STRUCTOPT

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor,
	   class ProblemDescriptor,
	   class ComplexVecType >
void
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optPreSolve()
{
  this->probDesc->preoptProcess();
  
  //this->allOps        = this->probDesc->getSolver();
  this->postProcessor = this->probDesc->getPostProcessor();

  int numLC = this->probDesc->getNumLC();
  
  solLC   = new VecType*[numLC];
  rhsLC   = new VecType*[numLC];
  gradLC  = new VecType*[numLC];
  
  for(int lc=0; lc<numLC; ++lc) 
    {
      solLC[lc]   = new VecType(this->probDesc->solVecInfo());
      rhsLC[lc]   = new VecType(this->probDesc->solVecInfo());
      gradLC[lc]  = 0;
    }

  gradLC[0]  = new VecType(this->probDesc->solVecInfo());
  this->adj  = 0;
  return;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
void
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optSolve()
{
  this->probDesc->reBuild(); 
  int numLC = this->probDesc->getNumLC(); // number of load cases 
  for(int lc=0; lc<numLC; ++lc) 
    {    
      filePrint(stderr," ... Processing loadcase    %7d ...\n",lc+1);      
      this->probDesc->setActiveLC(lc);      
      this->probDesc->getRHS(*rhsLC[lc]);

      this->probDesc->densProjectVector(*rhsLC[lc]);
      this->allOps->solve(*rhsLC[lc],*solLC[lc]);
      this->probDesc->densProjectVector(*solLC[lc]);
    }
  return;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
void
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optSolveDuds()
{
  int numLC = this->probDesc->getNumLC();  
  for(int lc=0; lc<numLC; ++lc) 
    {
      if (!gradLC[lc])  gradLC[lc]  = new VecType(this->probDesc->solVecInfo());       
      this->probDesc->setActiveLC(lc);
      this->probDesc->getPseudoLoad(*gradLC[lc],*solLC[lc]);      
      this->probDesc->densProjectVector(*gradLC[lc]);
      this->allOps->reSolve(*gradLC[lc]);
      this->probDesc->densProjectVector(*gradLC[lc]);
    }
  return;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
void
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optAdjSolve()
{
  adjNorm = (*this->adj).norm();
  if(adjNorm) 
    {
      this->probDesc->densProjectVector(*this->adj);
      this->allOps->reSolve(*this->adj);
      this->probDesc->densProjectVector(*this->adj);
    }
  return;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
double
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optGetAdjointPseudo(int lc)
{
  double val=0; 
  if (adjNorm) 
    { 
      val = this->probDesc->getAdjPseudoLoad(*solLC[lc],*this->adj);
    }  
 return val;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
VecType * 
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::getpsol (int lc) 
{
  // check if load case exists 
  int numLC = this->probDesc->getNumLC();  
  if ( lc >= numLC ) 
    {
      fprintf(stderr," *** ERROR: load case id %d > number of load cases %d\n",
	      lc+1,numLC);
      exit(-1);
    } 
  return solLC[lc];
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
VecType *  
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::getpgrad (int lc) 
{  
  // check if load case exists 
  int numLC = this->probDesc->getNumLC();  
  if(lc >= numLC) 
    {
      fprintf(stderr," *** ERROR: load case id %d > number of load cases %d\n",
	      lc+1,numLC);
      exit(-1);
    } 
  return gradLC[lc];
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
void
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::optPostProcessor(double time)
{
  int numLC = this->probDesc->getNumLC();  
  for(int lc=0; lc<numLC; ++lc) 
    {    
      this->probDesc->setActiveLC(lc);      
      this->postProcessor->staticOutput(*solLC[lc], *rhsLC[lc], time);
    }
  return;
}

//------------------------------------------------------------------------------
template < class Scalar,
	   class OpSolver, 
	   class VecType, 
	   class PostProcessor, 
	   class ProblemDescriptor,
	   class ComplexVecType >
VecType*
StaticSolver_opt< Scalar, OpSolver, VecType, 
		  PostProcessor, ProblemDescriptor, ComplexVecType >
::getpadj()
{
  if ( ! this->adj ) 
    {
      this->adj  = new VecType(this->probDesc->solVecInfo());  
    }
  return this->adj;
}

#endif //#ifdef STRUCTOPT
