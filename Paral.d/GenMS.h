#ifndef _GEN_M_S_H_
#define _GEN_M_S_H_

template<class Scalar>
class GenMultiSparse : public GenSparseMatrix<Scalar> 
{
   GenSparseMatrix<Scalar> *K, *Kii;
   GenAssembledFullM<Scalar> *Kcc;
   GenCuCSparse<Scalar>      *Krc;
   GenCuCSparse<Scalar>      *Kib;
   GenDBSparseMatrix<Scalar> *Kbb;
  public:
   GenMultiSparse(GenSparseMatrix<Scalar> *_K, GenSparseMatrix<Scalar> *_Kii, 
                  GenDBSparseMatrix<Scalar> *_Kbb, GenCuCSparse<Scalar> *_Kib) 
   { K= _K; Kii = _Kii; Kbb = _Kbb; Kib = _Kib; }

   GenMultiSparse(GenSparseMatrix<Scalar> *_K, GenSparseMatrix<Scalar> *_Kii, 
                  GenDBSparseMatrix<Scalar> *_Kbb, GenCuCSparse<Scalar> *_Kib,
                  GenCuCSparse<Scalar> *_Krc, GenAssembledFullM<Scalar> *_Kcc)
   { K= _K; Kii = _Kii; Kbb = _Kbb; Kib = _Kib; Krc = _Krc; Kcc = _Kcc; }
   
   ~GenMultiSparse() {};
    
   void add(FullSquareMatrix & kel, int *dofs);
   void addDiscreteMass(int dof, Scalar mass);
   Scalar diag(int i) const { return K->diag(i); }
   Scalar &diag(int i) { return K->diag(i); }
   int dim() { return K->dim() ; }
   int neqs() { return K->neqs() ; }
   void zeroAll();
};

template<class Scalar>
void
GenMultiSparse<Scalar>::zeroAll() 
{
 K->zeroAll();
 if(Kbb) Kbb->zeroAll();
 if(Kib) Kib->zeroAll();
 if(Kii) Kii->zeroAll();
 if(Kcc) Kcc->zero();
 if(Krc) Krc->zeroAll();
}

template<class Scalar>
void
GenMultiSparse<Scalar>::add(FullSquareMatrix & kel, int *dofs) 
{
 if(K)     K->add(kel, dofs);
 if(Krc) Krc->add(kel, dofs);
 if(Kcc) Kcc->add(kel, dofs);
 if(Kbb) Kbb->add(kel, dofs);
 if(Kib) Kib->add(kel, dofs);
 if(Kii) Kii->add(kel, dofs);
}

template<class Scalar>
void
GenMultiSparse<Scalar>::addDiscreteMass(int dof, Scalar mass)
{
// Check on if DOF is unconstrained must be done
// before this is called.
 if(K)     K->addDiscreteMass(dof, mass);
 if(Kbb) Kbb->addDiscreteMass(dof, mass);
 if(Kib) Kib->addDiscreteMass(dof, mass);
 if(Kii) Kii->addDiscreteMass(dof, mass);
 if(Kcc) Kcc->addDiscreteMass(dof, mass);
 if(Krc) Krc->addDiscreteMass(dof, mass);
}

//-----------------------------------------------------------------

#endif
