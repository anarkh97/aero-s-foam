template<class Scalar>
void
GenSubDomain<Scalar>::multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw) const
{
  // localvec += Br^T * interfvec
  int i, iDof, k;
  bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(i = 0; i < numMPC; ++i) mpcFlag[i] = true;

  for(iDof = 0; iDof < totalInterfSize; ++iDof) {
    switch(boundDofFlag[iDof]) {
      case 0:
        localvec[allBoundDofs[iDof]] += interfvec[iDof];
        break;
      case 1:  // active wet interface
        uw[-1-allBoundDofs[iDof]] = -interfvec[iDof];
        break;
      case 2: { // dual mpc
        int locMpcNb = -1-allBoundDofs[iDof];
        if(mpcFlag[locMpcNb]) { 
          SubLMPCons<Scalar> *m = mpc[locMpcNb];
          for(k = 0; k < m->nterms; k++) {
            int dof =(m->terms)[k].ccdof;
            Scalar coef = (m->terms)[k].coef;
            if(dof >= 0) localvec[dof] += interfvec[iDof]*coef;
          }
          mpcFlag[locMpcNb] = false;
        }
      } break;
    }
  }
 
 // coupled_dph: localvec -= Krw * uw
 if(Krw) Krw->multAddNew(uw, localvec);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multBr(const Scalar *localvec, Scalar *interfvec, Scalar *uc, Scalar *uw) const
{
 // interfvec = Br * localvec
 int iDof, k;
 for(iDof = 0; iDof < totalInterfSize; ++iDof) {
   switch(boundDofFlag[iDof]) {
     case 0:
       interfvec[iDof] = localvec[allBoundDofs[iDof]];
       break;
     case 1:   // active wet interface
       interfvec[iDof] = -uw[-1-allBoundDofs[iDof]];
       break;
     case 2: { // dual mpc
       int locMpcNb = -1-allBoundDofs[iDof];
       SubLMPCons<Scalar> *m = mpc[locMpcNb];
       interfvec[iDof] = 0;
       for(k = 0; k < m->nterms; k++) {
         int cc_dof = (m->terms)[k].ccdof;
         Scalar coef = (m->terms)[k].coef;
         if(cc_dof >= 0) interfvec[iDof] += localvec[cc_dof] * coef;  // not a corner
         else {
           int c_dof = (m->terms)[k].cdof;
           if(c_dof >= 0) { // corner or wet interface
             int dof = (m->terms)[k].dof;
             if((dof >= 0) && (cornerMap[dof] >= 0)) {  // corner
               if(cornerEqNums[cornerMap[dof]] >= 0)
                 interfvec[iDof] += uc[cornerEqNums[cornerMap[dof]]] * coef;
             }
           }
         }
       }
     } break;
   }
 }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multAddCT(const Scalar *interfvec, Scalar *localvec) const
{
  // localvec += C^T * interfvec
  bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
  for(int i=0; i<numMPC; ++i) mpcFlag[i] = true;

  for(int i = 0; i < scomm->lenT(SComm::mpc); i++) {
    int locMpcNb = scomm->mpcNb(i);
    if(mpcFlag[locMpcNb]) {
      SubLMPCons<Scalar> *m = mpc[locMpcNb];
      int iDof = scomm->mapT(SComm::mpc,i);
      for(int k = 0; k < m->nterms; k++) {
        int c_dof = (m->terms)[k].cdof;
        if(c_dof >= 0) localvec[c_dof] += interfvec[iDof]*(m->terms)[k].coef;
      }
      mpcFlag[locMpcNb] = false;
    }
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multC(const Scalar *localvec, Scalar *interfvec) const
{
  // interfvec = C * localvec
  for(int i = 0; i < scomm->lenT(SComm::mpc); i++) {
    int locMpcNb = scomm->mpcNb(i);
    SubLMPCons<Scalar> *m = mpc[locMpcNb];
    int iDof = scomm->mapT(SComm::mpc,i);
    interfvec[iDof] = 0;
    for(int k = 0; k < m->nterms; k++) {
      int c_dof = (m->terms)[k].cdof;
      if(c_dof >= 0) interfvec[iDof] += localvec[c_dof]*(m->terms)[k].coef;
    }
  }
}

