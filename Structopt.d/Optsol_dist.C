#if defined(STRUCTOPT) && defined(DISTRIBUTED)

#include <Structopt.d/Optsol_dist.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optvar.h>
#include <cassert>
#include <algorithm>
#include <boost/scoped_array.hpp>

const int Optalg_dist::root = 0;

//------------------------------------------------------------------------------
void Optalg_dist::bcast_action(MPI_Comm comm, Optalg_action& a)
{
  int ierr = MPI_Bcast(&a, 1, MPI_INT, root, comm);
  assert(ierr == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist::bcast_int(MPI_Comm comm, int& i)
{
  int ierr = MPI_Bcast(&i, 1, MPI_INT, root, comm);
  assert(ierr == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------
template<>
void Optalg_dist::bcast_array(MPI_Comm comm, int* ia, size_t size)
{
  int ierr = MPI_Bcast(ia, size, MPI_INT, root, comm);
  assert(ierr == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------
template<>
void Optalg_dist::bcast_array(MPI_Comm comm, double* ia, size_t size)
{
  int ierr = MPI_Bcast(ia, size, MPI_DOUBLE, root, comm);
  assert(ierr == MPI_SUCCESS);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist::bcast_func_params(MPI_Comm comm, int& iter,
				    double* vars, size_t length)
{
  bcast_int(comm, iter);
  bcast_array(comm, vars, length);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist::bcast_grad_params(MPI_Comm comm,
				    int* active, int& actsize,
				    int& funcGrad,				    
				    double* vars, size_t varsize)
{
  bcast_int(comm, actsize);
  if(actsize >0) { bcast_array(comm, active, actsize); }
  bcast_int(comm, funcGrad);
  if(funcGrad == 0)
    { bcast_array(comm, vars, varsize); }
  return;
}


//------------------------------------------------------------------------------
void Optalg_dist_masterWrapper::solve(Optsol *_optsol)
{
  Optalg_action a = Optalg_begin;
  bcast_action(comm, a);
  // do actual solve here -- this is just a stub
  assert(0);
  a = Optalg_end;
  bcast_action(comm, a);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist_masterWrapper::func(int iter)
{
  Optalg_action a = Optalg_func;
  bcast_action(comm, a);
  bcast_func_params(comm, iter, getOptsol()->var, getOptsol()->numvar);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist_masterWrapper::grad(int* active, int funcGrad)
{
  Optalg_action a = Optalg_grad;
  bcast_action(comm, a);
  int numactive = (active==0)? 0 : std::max(getOptsol()->numcon, 1);
  bcast_grad_params(comm, active, numactive, funcGrad, 
		    getOptsol()->var, getOptsol()->numvar);
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist_slaveFacade::solve(Optsol *_optsol)
{
  optsol = _optsol;
  
  Optalg_action a;
  int iter, size, funcGrad;
  boost::scoped_array<int> active(new int[std::max(getOptsol()->numcon, 1)]);
  
  bcast_action(comm, a);
  assert(a==Optalg_begin);
  while(a!=Optalg_end)
    {
      bcast_action(comm, a);
      switch(a)
	{
	case Optalg_end:
	  break;
	case Optalg_func:
	  bcast_func_params(comm, iter, 
			    getOptsol()->var, getOptsol()->numvar);
	  func(iter);
	  break;
	case Optalg_grad:
	  bcast_grad_params(comm, active.get(), size,
			    funcGrad, 
			    getOptsol()->var, getOptsol()->numvar);
	  if(size>0) 
	    { grad(active.get()); }
	  else 
	    { grad(); }
	  break;
	default:
	  assert(0);
	}
    }
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist_slaveFacade::func(int iter)
{
  getOptsol()->func(iter);   
  ordercon();  
  return;
}

//------------------------------------------------------------------------------
void Optalg_dist_slaveFacade::grad(int* active)
{
  int gtyp = optgrad.typ;
  reorderActive(active,getOptsol()->actcon);   
  if (!gtyp) { getOptsol()->optpro->optvar->updvar(); }
  optgrad.grad(getOptsol()->actcon);
  if (!gtyp) { getOptsol()->optpro->optvar->resetvar(); }
  ordergradcon();
  return;
}

//------------------------------------------------------------------------------
void OptalgNlpql_dist::solve(Optsol *_optsol)
{
  Optalg_action a = Optalg_begin;
  bcast_action(comm, a);
  OptalgNlpql::solve(_optsol);
  a = Optalg_end;
  bcast_action(comm, a);
  return;
}

//------------------------------------------------------------------------------
void OptalgNlpmma_dist::solve(Optsol *_optsol)
{
  Optalg_action a = Optalg_begin;
  bcast_action(comm, a);
  OptalgNlpmma::solve(_optsol);
  a = Optalg_end;
  bcast_action(comm, a);
  return;
}

/* PJSA
//------------------------------------------------------------------------------
void OptalgSnopt_dist::solve(Optsol *_optsol)
{
  Optalg_action a = Optalg_begin;
  bcast_action(comm, a);
  OptalgSnopt::solve(_optsol);
  a = Optalg_end;
  bcast_action(comm, a);
  return;
}
*/

#endif
