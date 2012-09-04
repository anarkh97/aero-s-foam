#include <Structopt.d/Structopt_dec.h>

//------------------------------------------------------------------------------
template<class Scalar>
void GenDistrDomain_opt<Scalar>::preProcess() 
{
  GenDistrDomain<Scalar>::preProcess();
  assert(this->subToNode != 0);  
  if(this->nodeToSub == 0)
    { this->nodeToSub = this->subToNode->reverse(); }
  assert(this->nodeToSub != 0);  
  for(int i=0; i<this->getNumSub(); ++i)
    {
      assert(dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i)) != 0);
      dynamic_cast<GenSubDomain_opt<Scalar>*>(this->getSubDomain(i))->makeGlobalToLocalElementMap(); 
    }
  return;
}

/*
//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getStructureMass(int* eleList, int listSize)
{
  double locMass = GenDecDomain_opt<Scalar>::getStructureMass(eleList, listSize);
  this->communicator->globalSum(1, &locMass);
  return locMass;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getStrainEnergy(GenDistrVector<Scalar>& sol, int* eleList, int listSize)
{
  double locEn = GenDecDomain_opt<Scalar>::getStrainEnergy(sol, eleList, listSize);
  this->communicator->globalSum(1, &locEn);
  return locEn;  
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getGradStrainEnergy(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad,
						       int* eleList, int listSize)
{
  double locdEn = GenDecDomain_opt<Scalar>::getGradStrainEnergy(sol, grad, eleList, listSize);
  this->communicator->globalSum(1, &locdEn);
  return locdEn;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getGradPartStrainEnergy(GenDistrVector<Scalar>& sol, 
							   int* eleList, int listSize)
{
  double locdEn = GenDecDomain_opt<Scalar>::getGradPartStrainEnergy(sol, eleList, listSize);
  this->communicator->globalSum(1, &locdEn);
  return locdEn;  
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getGradStructureMass(int* eleList, int listSize)
{
  double locdMass = GenDecDomain_opt<Scalar>::getGradStructureMass(eleList, listSize);
  this->communicator->globalSum(1, &locdMass);
  return locdMass;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDistrDomain_opt<Scalar>::DStiffnesDSmDUAdj(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& adj)
{
  Scalar result = GenDecDomain_opt<Scalar>::DStiffnesDSmDUAdj(sol, adj);
  this->communicator->globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
void GenDistrDomain_opt<Scalar>::getMidPointMass(double& totmas, double* midPoint)
{
  GenDecDomain_opt<Scalar>::getMidPointMass(totmas, midPoint);
  this->communicator->globalSum(1, &totmas);
  this->communicator->globalSum(3, midPoint);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDistrDomain_opt<Scalar>::getNodalDisp(GenDistrVector<Scalar>& sol, int node, int dispTyp)
{
  Scalar result = GenDecDomain_opt<Scalar>::getNodalDisp(sol, node, dispTyp);
  this->communicator->globalSum(1, &result);
  return result;
}


//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDistrDomain_opt<Scalar>::getGradNodalDisp(GenDistrVector<Scalar>& sol, GenDistrVector<Scalar>& grad, 
						    int node, int dispTyp)
{
  Scalar result = GenDecDomain_opt<Scalar>::getGradNodalDisp(sol, grad, node, dispTyp);
  this->communicator->globalSum(1, &result);
  return result;
}


//------------------------------------------------------------------------------
template<class Scalar>
Scalar GenDistrDomain_opt<Scalar>::getAverageDisp(GenDistrVector<Scalar>& sol,
						const int* nodes, const int* dtypes,
						int size)
{
  Scalar result = GenDecDomain_opt<Scalar>::getAverageDisp(sol, nodes, dtypes, size);
  this->communicator->globalSum(1, &result);
  return result;
}

//------------------------------------------------------------------------------
template<class Scalar>
double GenDistrDomain_opt<Scalar>::getDisplevel(GenDistrVector<Scalar>& sol, 
					      int aFlag, int quoFlag, int size, 
					      const int* nodeid, const int* distyp, 
					      const double* refVal, const double* difVal, double powFac)
{
  double result = pow(GenDecDomain_opt<Scalar>::getDisplevel(sol, aFlag, quoFlag, size, nodeid,
							     distyp, refVal, difVal, powFac),
		      powFac);
  this->communicator->globalSum(1, &result);
  return pow(result, 1.0/powFac);
}
*/
