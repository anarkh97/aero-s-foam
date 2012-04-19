#ifndef _MODAL_BASE_H_
#define _MODAL_BASE_H_

#include <Math.d/matrix.h>
#include <Driver.d/Domain.h>

class Domain;
class PrevFrc;
template <class V> class SysState;

class DiagonalMatrix
{
/* a minimalist class for a diagonal matrix
   used as the data type for the members of ModalOps
   numRbm() is present to satisfy call from DynamProbType::solve()
*/
private:

  int neq;    // ie number of equations, the number of diagonal entries
  double *d;  // the diagonal entries

public:

  DiagonalMatrix(){ d = 0; }
  DiagonalMatrix(int b){ neq = b; d = new double[b]; }
  ~DiagonalMatrix(){ if (d) delete [] d; d = 0; }

  double &operator[](int i) const { return d[i]; }

  void setDiag(double val);
  void mult(Vector &v, Vector &Av);
  void invertDiag();
  void reSolve(Vector &rhs);

  int numRBM(){ return 0; }
};

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

struct ModalOps
{
/* matrices for ModalDescr
   M     mass matrix
   C     damping matrix
   K     stiffness matrix
 dynMat  dynamic operator matrix for time integration
*/
  DiagonalMatrix *M, *C, *K, *dynMat, *Msolver;
  SparseMatrix *Muc, *Cuc;
};

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

class ModalBase
{
/* base class for NLModalDescr
   NOTE to self: if modalizing dsp/vel/acc is desired,
     will need to store scale as data memeber
*/
protected:

  Domain *domain;
  Vector *modesFl;   // eigenmodes of non-zero freq
  Vector *modesRB;   // rigid body modes
  double *freqs;     // circular frequencies (rad/s) of modesFl
  int numModes, numRBM, numFlex, numConstr;
  int *cDofIdx;      // stores index of the dof for each constraint

  double mass;    // total mass of the structure
  double cg[3];   // location of the center of gravity

  double *bcx, *vcx;    // prescribed displacement and velocity values

  Vector fullTmpF, fullTmpGrav, fullAeroF;
  Vector fullDsp, fullVel, fullAcc, fullPrevVel;
  PrevFrc *prevFrc, *prevFrcBackup;

public:

  virtual ~ModalBase() {/* TODO */};  
  ModalBase(){}
  ModalBase(Domain *d){ domain = d; }

  virtual void preProcess() = 0;
  void preProcessBase();
  void populateRBModes();
  void populateFlexModes(double scale = 1.0, bool readAll = 0);

  void initStateBase(Vector& dsp, Vector& vel, Vector& acc,
    Vector& vel_p, int idxOffset = 0);

  void outputModal(SysState<Vector> &state, Vector& extF, int tIndex);
};

#endif
