#include <Utils.d/dbg_alloca.h>
#include <stdio.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/Mpc.h>

SparseData::~SparseData()
{
  clean_up();
}

void
SparseData::clean_up()
{
  if(unconstrNum && myMem) { delete [] unconstrNum; unconstrNum=0; }
  if(constrndNum && myMem) { delete [] constrndNum; constrndNum=0; }
  if(xunonz)               { delete [] xunonz;      xunonz=0; }
  if(rowu && myMem_rowu)   { delete [] rowu;        rowu=0; }
  if(colu)		   { delete [] colu;	    colu=0; }	
}

SparseData::SparseData()
{
  initialize();
}

void
SparseData::initialize()
{
   unconstrNum    = 0;
   constrndNum    = 0;
   xunonz         = 0;
   rowu           = 0;
   colu		  = 0;
   numConstrained = 0;
   numUncon       = 0;
   neq            = 0;
   myMem          = 0;
   myMem_rowu     = 1;
}

SparseData::SparseData(Connectivity *con, DofSetArray *dsa, int *bc)
{
 initialize();
 int i, j, k, l;

 neq = dsa->size(); // Total number of equations
 int numNodes = dsa->numNodes(); // number of nodes
 if(con->csize() < numNodes) numNodes = con->csize();

 // We build a temporary dof to Node table.
 myMem = 1;
 unconstrNum = new int[neq];
 constrndNum = new int[neq];

 int cn = 0, un = 0;
 for(i=0; i<neq; ++i) {
   if(bc[i] == BCFIXED) {
     unconstrNum[i] = -1;
     constrndNum[i] = cn;
     cn = cn+1;
   }
   else {
     unconstrNum[i] = un;
     constrndNum[i] = -1;
     un = un+1;
   }
 }
 
 numConstrained = cn;

 int *numFixed    = (int * ) dbg_alloca(sizeof(int)*numNodes);
 int *numNotFixed = (int * ) dbg_alloca(sizeof(int)*numNodes);

 for(i=0; i<numNodes; ++i) {
   numFixed[i] = 0; 
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j=0; j<myNumDofs; ++j)
     if(constrndNum[myFirstDof+j] >= 0) numFixed[i] += 1;
   numNotFixed[i] = myNumDofs - numFixed[i];  
 }

 // Allocate memory for xunonz
 xunonz = new int[numConstrained+1];
 
 for(i=0; i<numNodes; ++i) {
   if(numFixed[i] == 0) continue;
   int nentries=0;
   for(j=0; j < con->num(i); ++j) 
     nentries += numNotFixed[(*con)[i][j]];
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j=0; j<myNumDofs; ++j)
     if(constrndNum[myFirstDof+j] >= 0)
        xunonz[constrndNum[myFirstDof+j]] = nentries;
 }

 int count = 0;
 for(i=0; i<numConstrained; ++i) {
   int temp = xunonz[i];
   xunonz[i] = count;
   count += temp;
 }
 // count = # of entries in kuc

 xunonz[numConstrained] = count;

// Allocate memory for rowu i.e. row #s

 rowu = new int[count];

 for(i=0; i<numNodes; ++i) {
   if(numFixed[i] == 0) continue;
   int iFirstDof = dsa->firstdof(i);
   int nentries = 0;
   for(j =0; j < con->num(i); ++j) {
     int jNode = (*con)[i][j];
     if(numNotFixed[jNode] == 0) continue;
     int jFirstDof = dsa->firstdof(jNode);
     for(k = 0; k < dsa->weight(jNode); ++k) {
       if(unconstrNum[jFirstDof + k] >= 0) {
         for(l = 0; l < dsa->weight(i); ++l)
            if(constrndNum[iFirstDof+l] >= 0)
               rowu[ xunonz [constrndNum[iFirstDof+l]] +nentries] = unconstrNum[jFirstDof + k];
         nentries++; 
       }
     }
   }
 }

}

SparseData::SparseData(Connectivity *con, DofSetArray *dsa, 
                       DofSetArray *c_dsa)
{
 initialize();
 int i, j, k, l;

 neq = dsa->size();

 int numNodes = c_dsa->numNodes();
 if(con->csize() < numNodes) numNodes = con->csize();

 numUncon       = c_dsa->size();  // number of unconstrained dof
 myMem = 0;
 unconstrNum = c_dsa->getUnconstrNum();
 if(numUncon == 0) return;
 numConstrained = neq - numUncon; // number of constrained dof
 constrndNum = c_dsa->getConstrndNum();
 if(numConstrained == 0) return;
 
 int *numFixed    = (int * ) dbg_alloca(sizeof(int)*numNodes);
 int *numNotFixed = (int * ) dbg_alloca(sizeof(int)*numNodes);

 for(i=0; i<numNodes; ++i) {
   numFixed[i] = 0;
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j = 0; j < myNumDofs; ++j)
     if(constrndNum[myFirstDof+j] >= 0) numFixed[i] += 1;
   numNotFixed[i] = myNumDofs - numFixed[i];
 }

// Allocate memory for xunonz.
 xunonz = new int[numConstrained+1];

 for(i=0; i<numNodes; ++i) {
   if(numFixed[i] == 0) continue;
   int nentries=0;
   for(j=0; j<con->num(i); ++j)
     nentries += numNotFixed[(*con)[i][j]];
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j=0; j<myNumDofs; ++j)
     if(constrndNum[myFirstDof+j] >= 0)
        xunonz[constrndNum[myFirstDof+j]] = nentries;
 }

 int count = 0;
 for(i=0; i<numConstrained; ++i) {
   int temp = xunonz[i];
   xunonz[i] = count;
   count += temp;
 }
 xunonz[numConstrained] = count;

 rowu = new int[count];

 for(i=0; i<numNodes; ++i) {
   if(numFixed[i] == 0) continue;
   int iFirstDof = dsa->firstdof(i);
   int nentries=0;
   for(j=0; j<con->num(i); ++j) {
     int jNode = (*con)[i][j];
     if(numNotFixed[jNode] == 0) continue;
     int jFirstDof = dsa->firstdof(jNode);
     for(k=0; k<dsa->weight(jNode); ++k) {
       if(unconstrNum[jFirstDof + k] >= 0) {
         for(l = 0; l < dsa->weight(i); ++l)
            if(constrndNum[iFirstDof+l] >= 0)
               rowu[ xunonz [constrndNum[iFirstDof+l]] +nentries] = unconstrNum[jFirstDof + k];
         nentries++;
       }
     }
   }
 }

}

SparseData::SparseData(Connectivity *con, DofSetArray *dsa, int *glBoundMap, int *glInternalMap)
{
 initialize();
 int i,j,k,l,iDof;

 // Get the global length
 neq = dsa->size();
 myMem = 0;
 constrndNum = glBoundMap;
 unconstrNum = glInternalMap;
 // Get the number of nodes
 int numNodes = dsa->numNodes();
 if(con->csize() < numNodes) numNodes = con->csize();

 // Compute number of Boundary dofs and number of Internal dofs
 int boundaryLen = 0;
 int internalLen = 0;
 for(iDof = 0; iDof < neq; ++iDof) {
   if(glBoundMap[iDof] >= 0 ) boundaryLen += 1;
   if(glInternalMap[iDof] >= 0) internalLen += 1;
 }

// Allocate memory for xunonz
 xunonz  = new int[boundaryLen+1];
 numConstrained = boundaryLen;
 numUncon       = internalLen; //HB
 int *numBound    = (int * ) dbg_alloca(sizeof(int)*numNodes);
 int *numInternal = (int * ) dbg_alloca(sizeof(int)*numNodes);

 for(i=0; i<numNodes; ++i) {
   numBound[i]   = 0;
   numInternal[i] = 0;
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j = 0; j < myNumDofs; ++j){
     if(glBoundMap[myFirstDof+j] >= 0) numBound[i] += 1;
     if(glInternalMap[myFirstDof+j] >= 0) numInternal[i] += 1;
   }
 }

 for(i=0; i<numNodes; ++i) {
   if(numBound[i] == 0) continue;
   int nentries=0;
   for(j=0; j<con->num(i); ++j)
     nentries += numInternal[(*con)[i][j]];
   int myFirstDof = dsa->firstdof(i);
   int myNumDofs  = dsa->weight(i);
   for(j=0; j<myNumDofs; ++j)
     if(glBoundMap[myFirstDof+j] >= 0)
        xunonz[glBoundMap[myFirstDof+j]] = nentries;
 }

 int count = 0;
 for(i=0; i<boundaryLen; ++i) {
   int temp = xunonz[i];
   xunonz[i] = count;
   count += temp;
 }
 xunonz[boundaryLen] = count;

 rowu = new int[count];
 for(i=0; i<numNodes; ++i) {
   if(numBound[i] == 0) continue;
   int iFirstDof = dsa->firstdof(i);
   int nentries=0;
   for(j=0; j<con->num(i); ++j) {
     int jNode = (*con)[i][j];
     if(numInternal[jNode] == 0) continue;
     int jFirstDof = dsa->firstdof(jNode);
     for(k=0; k<dsa->weight(jNode); ++k) {
       if(glInternalMap[jFirstDof + k] >= 0) {
         for(l = 0; l < dsa->weight(i); ++l)
            if(glBoundMap[iFirstDof+l] >= 0)
               rowu[ xunonz [glBoundMap[iFirstDof+l]] +nentries]
                   = glInternalMap[jFirstDof + k];
         nentries++;
       }
     }
   }
 }
}

// used for Mumps & Spooles (1st constructor)
// Constructor for DBSparseMatrix data structures
SparseData::SparseData(EqNumberer *_dsa, Connectivity *cn, int* rCN, int expand, int make_colu)
{
  initialize();
  int i, j, k, thisNode;

  neq = _dsa->size();

// We build a temporary dof to Node table.
  int* dofToN = new int[neq]; //HB
  int numNodes = _dsa->numNodes();

  int *firstDOF   = (int *) dbg_alloca(sizeof(int)*numNodes);
  int *nodeWeight = (int *) dbg_alloca(sizeof(int)*numNodes);

  numUncon  = 0;
  myMem = 1;
  unconstrNum = new int[neq];
  for(i = 0; i < neq; ++i) {
    unconstrNum[i] = (rCN) ? rCN[i] : i;
    if(unconstrNum[i] >= 0) numUncon++;
  }
  //cerr << "cn = \n"; cn->print(); "rCN = "; for(int i=0;i<neq;++i) cerr << rCN[i] << " "; cerr << endl;

  for(i=0; i < numNodes; ++i) {
    int fdof = _dsa->firstdof(i);
    int ndof = _dsa->weight(i);
    nodeWeight[i] = 0;
    for(j=0; j<ndof; ++j) {
      dofToN[j+fdof] = i;
      if(unconstrNum[j+fdof] >= 0) {
             nodeWeight[i]++;
      }
    }
  }

  for(i=0; i < numNodes; ++i)
    firstDOF[i] = -1;

  for(i = 0; i < neq; ++i) {
    int lDof = unconstrNum[i];
    if(lDof >= 0) {
      int lNode = dofToN[i];
      if(firstDOF[lNode] < 0 || lDof < firstDOF[lNode])
        firstDOF[lNode]  = lDof;
    }
  }

  // We assume that any Dof that has been used has a diagonal term.
  xunonz = new int[numUncon+1];
  xunonz[0] = 1;

  for(i=0; i<neq; ++i) {
     k = unconstrNum[i];     //  or k = unconstrndNum[i]
     if(k == -1) continue;
     xunonz[k+1]  = xunonz[k];
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
       int jnode = (*cn)[thisNode][j];
       int fjdof = firstDOF[jnode];
       if(fjdof < 0 || fjdof > k) continue;
       if(nodeWeight[jnode] + fjdof <= k)
         xunonz[k+1] += nodeWeight[jnode];
       else
         xunonz[k+1] += k - fjdof + 1;
     }
  }

// Allocate memory for rowu (row numbers)
  rowu  = new int[xunonz[numUncon]];
  for(i=0; i<neq; ++i) {
     int m = unconstrNum[i];
     if(m == -1) continue;
     int numFound = 0;
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
       int jnode = (*cn)[thisNode][j];
       int fjdof = firstDOF[jnode];
       if(fjdof < 0 || fjdof > m) continue;
       for(k=0; k<nodeWeight[jnode]; ++k) 
         if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
     }
     rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.  
  }

/* KAS.
   Construct colu : column numbers for each element - used for Mumps
   same size as rowu. use xunonz index to get column number.*/
  if(make_colu) { 
     colu  = new int[xunonz[numUncon]];
     k     = 0;
     for (i = 0; i < numUncon; i++) {
        // number of off-diagonal elements, will repeat column index numOffDiag times in colu
        int numOffDiag   = xunonz[i+1] - xunonz[i];   // xunonz is of size numUncon +1 
	for (j = 0; j < numOffDiag; j++) {
           colu[k++]  = i+1;
        }
     }
  }
  if(dofToN) delete [] dofToN; //HB
}

// used by 2nd Mumps constructor (with make_colu = 1) and 2nd Spooles constructor
SparseData::SparseData(DofSetArray *_dsa, DofSetArray *c_dsa,
                       Connectivity *cn, int expand, int make_colu)
{
  initialize();
  int i, j, k, thisNode;

  neq = _dsa->size(); // Total number of equations
  numUncon = c_dsa->size(); // Number of unconstrained dof (equations)

  // We build a temporary dof to Node table.
  int* dofToN = new int[neq]; //HB
  int numNodes = _dsa->numNodes();

  for(i=0; i < numNodes; ++i) {
    int fdof = _dsa->firstdof(i);
    int ndof = _dsa->weight(i);
    for(j=0; j<ndof; ++j)
      dofToN[j+fdof] = i;
  }
  
  myMem = 0;
  unconstrNum = c_dsa->getUnconstrNum();

  // We assume that any Dof that has been used has a diagonal term.
  // Allocate memory for the diagonal location pointers
  xunonz = new int[numUncon+1];

  // Set the diagonal location pointers
  xunonz[0] = 1;
  for(i=0; i<neq; ++i) {
    k = unconstrNum[i]; 
    if(k == -1) continue;
    xunonz[k+1] = xunonz[k];
    thisNode = dofToN[i];
    for(j=0; j<cn->num(thisNode); ++j) {
      int jnode = (*cn)[thisNode][j];
      int fjdof = c_dsa->firstdof(jnode);
      if(fjdof < 0 || fjdof > k) continue;
      if(c_dsa->weight(jnode) + fjdof <= k)
        xunonz[k+1] += c_dsa->weight(jnode);
      else
        xunonz[k+1] += k - fjdof + 1;
    }
  }

  // Allocate memory for rowu (row numbers)
  rowu  = new int[xunonz[numUncon]];

  // Set the row numbers
  for(i=0; i<neq; ++i) {
    int m = unconstrNum[i];
    if(m == -1) continue;
    int numFound = 0;
    thisNode = dofToN[i];
    for(j=0; j<cn->num(thisNode); ++j) {
      int jnode = (*cn)[thisNode][j];
      int fjdof = c_dsa->firstdof(jnode);
      if(fjdof < 0 || fjdof > m) continue;
      for(k=0; k<c_dsa->weight(jnode); ++k)
        if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
    }
    rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.
  }
  
  if(dofToN) delete [] dofToN; //HB

  if(expand && numUncon > 0) { // PJSA
    int k,j;
    int *new_xunonz = new int[numUncon+1];
    int* counter = new int[numUncon+1]; //HB
    for(k=0; k < numUncon+1; k++) {
      new_xunonz[k] = xunonz[k];
      counter[k] = 0;
    }
    for(k=0; k < numUncon; k++)
      for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
        if(rowu[j]-1 != k)
          counter[rowu[j]]++;
      }
    for(k=1; k < numUncon; k++)
      counter[k] += counter[k-1];
    for(k=1; k < numUncon; k++)
      new_xunonz[k] += counter[k];
    new_xunonz[numUncon] += counter[numUncon-1];

    for(k=0; k < numUncon; k++)
       counter[k] = 0;

    int *new_rowu  = new int[new_xunonz[numUncon]];

//  map the other side to rowu and unonz
    for(k=0; k < numUncon; k++) {
       for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
          new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
          counter[k]++;
          if(rowu[j]-1 != k) {
            new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
            counter[rowu[j]-1]++;
          }
       }
    }
   
    delete [] rowu;   rowu   = new_rowu;
    delete [] xunonz; xunonz = new_xunonz;
    delete [] counter; //HB
  }


/* KAS.
   Construct colu : column numbers for each element - used for Mumps
   same size as rowu. use xunonz index to get column number.*/
  if(make_colu) {
     colu  = new int[xunonz[numUncon]];
     k     = 0;
     for (i = 0; i < numUncon; i++) {
        // number of off-diagonal elements, will repeat column index numOffDiag times in colu
        int numOffDiag   = xunonz[i+1] - xunonz[i];   // xunonz is of siez numUncon +1
        for (j = 0; j < numOffDiag; j++) {
            colu[k++]  = i+1;
        }
     }
  }
}


SparseData::SparseData(DofSetArray *_dsa, int *glInternalMap,
                       Connectivity *cn, int expand)
{
  initialize();
  int i, j, k, thisNode;

  neq  =  _dsa->size(); // Total number of equations

  // Compute number of Boundary dofs and number of Internal dofs
  int internalLen = 0;
  int iDof;
  for(iDof = 0; iDof < neq; ++iDof) {
    if(glInternalMap[iDof] >= 0) internalLen += 1;
  }
  numUncon = internalLen;

  int numNodes = _dsa->numNodes();
  int *firstDOF   = (int *) dbg_alloca(sizeof(int)*numNodes);
  int *nodeWeight = (int *) dbg_alloca(sizeof(int)*numNodes);

// We build a temporary dof to Node table.
  int* dofToN = new int[neq]; //HB

  myMem=0;
  unconstrNum = glInternalMap;

  for(i=0; i < numNodes; ++i) {
    int fdof = _dsa->firstdof(i);
    int ndof = _dsa->weight(i);
    nodeWeight[i] = 0;
    for(j=0; j<ndof; ++j) {
      dofToN[j+fdof] = i;
      if(unconstrNum[j+fdof] >= 0) {
             nodeWeight[i]++;
      }
    }
  }

  for(i=0; i < numNodes; ++i)
    firstDOF[i] = -1;

  for(i = 0; i < neq; ++i) {
    int lDof = unconstrNum[i];
    if(lDof >= 0) {
      int lNode = dofToN[i];
      if(firstDOF[lNode] < 0 || lDof < firstDOF[lNode])
        firstDOF[lNode]  = lDof;
    }
  }

  // We assume that any Dof that has been used has a diagonal term.
  // Allocate memory for the diagonal location pointers
  xunonz = new int[numUncon+1];

  // Set the diagonal location pointers
  xunonz[0] = 1;
  for(i=0; i<neq; ++i) {
     k = unconstrNum[i]; 
     if(k == -1) continue;
     xunonz[k+1] = xunonz[k];
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
        int jnode = (*cn)[thisNode][j];
        int fjdof = firstDOF[jnode];
        if(fjdof < 0 || fjdof > k) continue;
        if(nodeWeight[jnode] + fjdof <= k)
           xunonz[k+1] += nodeWeight[jnode];
        else
           xunonz[k+1] += k - fjdof + 1;
     }
  }

  // Allocate memory for rowu (row numbers)
  rowu  = new int[xunonz[numUncon]];

  // Set the row numbers
  for(i=0; i<neq; ++i) {
     int m = unconstrNum[i];
     if(m == -1) continue;
     int numFound = 0;
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
        int jnode = (*cn)[thisNode][j];
        int fjdof = firstDOF[jnode];
        if(fjdof < 0 || fjdof > m) continue;
        for(k=0; k<nodeWeight[jnode]; ++k) {
          if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
        }
     }
     rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.
  }
  if(dofToN) delete [] dofToN; //HB

  if(expand) {
    int k,j;
    int *new_xunonz = new int[numUncon+1];
    int* counter = new int[numUncon+1]; //HB

    for(k=0; k < numUncon+1; k++) {
      new_xunonz[k] = xunonz[k];
      counter[k] = 0;
    }

    for(k=0; k < numUncon; k++)
      for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
        if(rowu[j]-1 != k)
          counter[rowu[j]]++;
      }

    for(k=1; k < numUncon; k++)
      counter[k] += counter[k-1];

    for(k=1; k < numUncon; k++)
      new_xunonz[k] += counter[k];

    new_xunonz[numUncon] += counter[numUncon-1];

    for(k=0; k < numUncon; k++)
      counter[k] = 0;

    int *new_rowu  = new int[new_xunonz[numUncon]];

    // map the other side to rowu and unonz
    for (k=0; k < numUncon; k++) {
       for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
          new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
          counter[k]++;
          if(rowu[j]-1 != k) {
            new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
            counter[rowu[j]-1]++;
          }
       }
    }

    delete [] rowu;   rowu   = new_rowu;
    delete []xunonz; xunonz = new_xunonz;
    delete [] counter;
  }

}

//#include <algo.h>
#include <algorithm> // PJSA: for sgi intel

SparseData::SparseData(EqNumberer *eqn, Connectivity *cn, double trbm)
{
 initialize();
 int i, j, k, l;

 numUncon = neq  = eqn->size(); // Total number of equations

 int numNodes = cn->csize();
 int *nodeCWeight = (int *) dbg_alloca(sizeof(int)*numNodes);
 
 for(i = 0; i < numNodes; ++i)
   // sort((*cn)[i], (*cn)[i+1]);
   std::sort((*cn)[i], (*cn)[i+1]);  // PJSA: for sgi intel

 for(i = 0; i < numNodes; ++i) {
    nodeCWeight[i] = 0;
    for(j = 0; j < cn->num(i); ++j)
      nodeCWeight[i] += eqn->weight((*cn)[i][j]);
 }

 // Allocate the column pointer 
 xunonz = new int[neq+1];
 
 for(i = 0; i < numNodes; ++i) {
   int fDof = eqn->firstdof(i);
   int nDof = eqn->weight(i);
   for(j = 0; j < nDof; ++j)
     xunonz[fDof+j] = nodeCWeight[i];
 }

 int count = 0;
 for(i = 0; i < neq; ++i) {
   int tmpc = count;
   count += xunonz[i];
   xunonz[i] = tmpc;
 } 
 xunonz[neq] = count;
 
 rowu = new int[count];

 for(i = 0; i < numNodes; ++i) {
   int fDof = eqn->firstdof(i);
   int nDof = eqn->weight(i);
   int index = 0;
   for(k = 0; k < cn->num(i); ++k) {
     int nodeNum = (*cn)[i][k];
     int fkDof = eqn->firstdof(nodeNum);
     int nkDof = eqn->weight(nodeNum);
     for(l = 0; l < nkDof; ++l, ++index)
        for(j = 0; j < nDof; ++j)
           rowu[ xunonz[fDof+j] + index] = fkDof+l;
   }
 }
}

// ----- EXPERIMENTAL -----------------------------------------------

SparseData::SparseData(Connectivity *cn, EqNumberer *eqn, double trbm,
                       int expand)
{
  initialize();
  int i, j, k, thisNode;

  neq  = eqn->size(); // Total number of equations
  numUncon = eqn->size();

  // Addition from UH 3-19-00
  myMem = 1;
  unconstrNum = new int[neq];
  constrndNum = new int[neq];

  int un = 0;
  for(i=0; i<neq; ++i) {
     unconstrNum[i] = un;
     constrndNum[i] = -1;
     un = un+1;
  }

  // We build a temporary dof to Node table.
  int numNodes = eqn->numNodes();
  if( cn->csize() < numNodes ) numNodes = cn->csize();

  int* dofToN = new int[neq]; //HB

  for(i=0; i < numNodes; ++i) {
    int fdof = eqn->firstdof(i);
    int ndof = eqn->weight(i);
    for(j=0; j<ndof; ++j)
      dofToN[j+fdof] = i;
  }

  // We assume that any Dof that has been used has a diagonal term.
  // Allocate memory for the diagonal location pointers
  xunonz = new int[neq+1];

  // Set the diagonal location pointers
  xunonz[0] = 1;

  for(i=0; i<neq; ++i) {
     xunonz[i+1] = xunonz[i];
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
        int jnode = (*cn)[thisNode][j];
        int fjdof = eqn->firstdof(jnode);
        if(fjdof < 0 || fjdof > i) continue;
        if( eqn->weight(jnode) + fjdof <= i )
           xunonz[i+1] += eqn->weight(jnode);
        else
           xunonz[i+1] += i - fjdof + 1;
     }
  }

  // Allocate memory for rowu (row numbers)
  rowu  = new int[xunonz[neq]];

  // Set the row numbers
  for(i=0; i<neq; ++i) {
     int numFound = 0;
     thisNode = dofToN[i];
     for(j=0; j<cn->num(thisNode); ++j) {
        int jnode = (*cn)[thisNode][j];
        int fjdof = eqn->firstdof(jnode);

        if(fjdof < 0 || fjdof > i) continue;

        for(k=0; k<eqn->weight(jnode); ++k)
          if(fjdof + k < i) rowu[xunonz[i]-1+numFound++] = fjdof + k + 1;
     }
     rowu[xunonz[i]-1+numFound++] = i + 1; // This is for the diagonal terms.
  }
  if(dofToN) delete [] dofToN;

  if(expand) {
    int k,j;
    int *new_xunonz = new int[numUncon+1];
    int* counter = new int[numUncon+1]; //HB

    for(k=0; k < numUncon+1; k++)  {
      new_xunonz[k] = xunonz[k];
      counter[k] = 0;
    }

    for(k=0; k < numUncon; k++)
       for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
           if(rowu[j]-1 != k)
              counter[rowu[j]]++;
       }

    for(k=1; k < numUncon; k++)
        counter[k] += counter[k-1];

    for(k=1; k < numUncon; k++)
        new_xunonz[k] += counter[k];

    new_xunonz[numUncon] += counter[numUncon-1];

    for(k=0; k < numUncon; k++)
       counter[k] = 0;

    int *new_rowu  = new int[new_xunonz[numUncon]];

    // map the other side to rowu and unonz
    for(k=0; k < numUncon; k++) {
       for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
           new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
           counter[k]++;

           if(rowu[j]-1 != k) {
              new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
              counter[rowu[j]-1]++;
           }
       }
    }

    delete [] rowu;   rowu   = new_rowu;
    delete [] xunonz; xunonz = new_xunonz;
    delete [] counter;
  }

  dbg_alloca(0);

}


SparseData::SparseData(LMPCons **mpc, int numMPC, DofSetArray *c_dsa)
{
  initialize();
  numConstrained = numMPC;
  numUncon = c_dsa->size();
  xunonz = new int[numConstrained+1];

  int i;
  xunonz[0] = 0;
  for(i=1; i<numConstrained+1; ++i) {
    xunonz[i] = xunonz[i-1] + mpc[i-1]->nterms;
  }

  rowu = new int[xunonz[numConstrained]];

  int mstart, mstop, m;
  for(i=0; i<numConstrained; ++i) {
    mstart = xunonz[i];
    mstop  = xunonz[i+1];
    for(m=mstart; m<mstop; ++m) {
      rowu[m] = c_dsa->locate(mpc[i]->terms[m].nnum,
                             (1 << mpc[i]->terms[m].dofnum));
    }
  }
}


SparseData::SparseData(int num, int *xyzCount, int *xyzList)
{
  initialize();
  numConstrained = num;
  xunonz = new int[numConstrained+1];

  xunonz[0] = 0;
  int i;
  for(i=0; i<numConstrained; ++i)
    xunonz[i+1] = xunonz[i] + xyzCount[i];

  myMem_rowu = 1;
  rowu = xyzList;
}

//KHP:to store local G as a rectangular sparse matrix

SparseData::SparseData(int numInterface,
                       int *glbmap, int numModes, int ldm)
{
  initialize();
  numConstrained = numModes;
  xunonz = new int[numConstrained+1];
  int i;
  xunonz[0] = 0;
  for(i=1; i<numConstrained+1; ++i) {
    xunonz[i] = xunonz[i-1] + numInterface;
  }
  int count = numInterface*numConstrained;
  rowu      = new int[count];

  int mstart, mstop, m;
  for(i=0; i<numConstrained; ++i) {
    mstart = xunonz[i];
    mstop  = xunonz[i+1];
    for(m=mstart; m<mstop; ++m) {
      rowu[m] = glbmap[m-mstart];
    }
  }
}

