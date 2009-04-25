#ifndef _SANDDOMAIN_H_
#define _SANDDOMAIN_H_

#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Paral.d/MDDynam.h>
//#include <Dist.d/DistDom.h>
#include <Interface.d/SandiaElem.h>
class MpcLocal;
enum face_type;  

template<class Scalar>
struct GenMatPair 
{
  GenSolver<Scalar> *smat;
  GenSparseMatrix<Scalar> *dmat;
  GenMatPair(GenSolver<Scalar> *_smat, GenSparseMatrix<Scalar> *_dmat)
      { smat = _smat; dmat = _dmat; }
};

template<class Scalar>
class GenSandiaSubD : public GenSubDomain<Scalar> 
{
  int *sandiaToLocalMap;
  Rbm *rbm;
  bool subSysMatrix;
  const int *map;
  int *subMap;
  unsigned short int *dofMap;

  long memoryK;    // subdomain stiffness memory
  long memoryPrec; // subdomain preconditioner memory
  int glNumMPC;

  int numContactNeighb;
  int *contactNeighb, *mpcNeighb;
  int *neighbPtr;
  Connectivity *ctcInterfaceDOFs, *mpcInterfaceDOFs;

 public:
  GenSandiaSubD(int _globalSubNumber, int _localSubNumber, int nnodes, const double *x, const double *y,
                const double *z, int *locToGlob, unsigned short int *dofMap, int numEle, int *eptr, int *econ, 
                const int *map, int *subMap = 0);
  virtual ~GenSandiaSubD();
  void initialize();
  void initSysMatrix(int neq, GenSolver<Scalar> *&smat, GenSparseMatrix<Scalar> *&diagmat);
  void zeroSysMatrix(GenSolver<Scalar> *smat, GenSparseMatrix<Scalar> *diagmat);
  void deleteSysMatrix(GenSolver<Scalar> *smat = 0, GenSparseMatrix<Scalar> *diagmat = 0);
  void addSysMatrix(int neq, const double *Kaa, const int *Kaj, const int *Kai, Scalar multiplier,
                    GenSolver<Scalar> *&smat, GenSparseMatrix<Scalar> *&diagmat, bool isK = true);
  void getRBMs(int neq, int nRBM, Scalar *salinasRBMs, Scalar *fetiRBMs);
  void makeF(int neq, Scalar *f, Scalar *lf);
  void getD(int neq, Scalar *d, Scalar *ld, bool assemble);
  void remap(int neq, Scalar *salinasDis, Scalar *fetiDis);
  double salinasToFetiVector(double *salinasDisp, double *fetiDis);
  double fetiToSalinasVector(double* fetiDis, double* salinasDis);
  Rbm *getRbm() { return rbm; }
  void setMpcGlobalTermIds(MpcLocal *mpclocal, int *global_node_nums, bool oneSubPerCPU);
  long getMemoryK()    { return memoryK;    }
  long getMemoryPrec() { return memoryPrec; }
  int Neq() { return this->c_dsa->size(); }
};

class FetiParams;

template<class Scalar>
class GenSandiaDomain: public GenDecDomain<Scalar> 
{
  const double *x, *y, *z;
  int num_nodes, num_elems, neq, ndofs_deleted, numads;
  int *elem_ptr, *elem_conn, *bc_node_dof, *procad, *com_ptr, *com_lis, *global_node_nums;
  const int *map, *bc_node_ids;
  unsigned short *dofmap_on_node;

  bool oneSubPerCPU;
  FetiInfo *finfo;
  int numberMpc;  // this is the global number of mpcs
  int numberCtc;
  GenSolver<Scalar> **sysMatrix;
  GenSparseMatrix<Scalar> **sparseK;
  MpcLocal *mpclocal;
  bool computeRbms;
  GenFetiDPSolver<Scalar> * fetiSolver;
  int *local_node_nums;
  int global_node_max;
  int mortar_type;
  double *initial_disp;
  //double *dx, *dy, *dz;
  int problem_type; 
  bool have_contact;
 
 public:
  GenSandiaDomain(Domain *d, const double *x, const double *y, const double *z,
                  int num_nodes, int *elem_ptr, int *elem_conn, int num_elems,
                  int neq, const int *map, int ndofs_deleted, const int *bc_node_ids,
                  int *bc_node_dof, int numads, int *procad, int *com_ptr, int *com_lis,
                  unsigned short *dofmap_on_node, int *global_node_nums, FetiParams *params);
  ~GenSandiaDomain();
  void makeSubD();
  void makeSComm();
  void constrain();
  void setMPCs(int numberMpc, MpcLocal *mpc);
  int  updateLocalMpcs(MpcLocal **mpc);
  void setWaveNumbers(FetiParams *params);
  void preProcess();
  void initSysMatrix();
  void zeroSysMatrix();
  void addSysMatrix(const double *Kaa, const int *Kaj,
                    const int *Kai, Scalar multiplier, bool isK = true);
  GenFetiDPSolver<Scalar> *getSolver();
  void updateSysMatrixInSolver();
  void makeF(Scalar *f, GenDistrVector<Scalar> *df);
  void getD(Scalar *d, GenDistrVector<Scalar> *dd, bool assemble = false);
  void getRBMs(Scalar *);
  int getNumberMpc() { return numberMpc; } // number of mpcs (global)
  void constraintProduct(int num_vect, const double* R[], Scalar** V, int trans);
  MpcLocal* getLocalMpcs() { return mpclocal; }
  void mpcForces(double *lambdas);
  void setSurfaceInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                              enum face_type *ftype, int *face_ptr, int *face_conn,
                              double (*face_normal)[3], double (*gap_vector)[3], int *glbl_surface_id,
                              double *normal_tol = 0, double *tangent_tol = 0);
  int setContactNormal(const double (*normal)[3]); // set the contact normal (i.e. usually n12[n])
  //int setContactGap(const double (*gapVectors)[3]); // set the normal gap (i.e. usually g12[n] = G12[n].n12[n]))
  int setContactGap(const double *disp = 0); 
  int updateContactGap(const double *DU); // g12[n,j] <- g12[n] - (Du2[n,j]-Du1[n,j]).n12[n]
                                          // where n = (time/load) step counter & j = (Newton) iteration counter
                                          // within the current (time/load) step n
  int getContactForces(double *ctcForces); // compute the (total) contact forces since the last call
                                           // to setContactSurfaces/setContactGap
  int numWaveDir() { return (finfo->dph_flag) ? finfo->numdir : 0; }

 private:
  void initialize();
  void setFetiInfo(FetiParams *params);
  void setSolInfo(FetiParams *params);
  void makeCpuToSub();
  void makeSubToSubEtc();
  void initSubSysMatrix(int iSub);
  void zeroSubSysMatrix(int iSub);
  void addSubSysMatrix(int neq, const double *Kaa, const int *Kaj,
                       const int *Kai, Scalar multiplier, bool isK);
  void makeSubF(int iSub, Scalar *f, GenDistrVector<Scalar> *df);
  void getSubD(int iSub, Scalar *d, GenDistrVector<Scalar> *dd, bool assemble);
  void subMpcForces(int iSub, double *lambdas);
  GenSandiaSubD<Scalar> * SubDomain(int iSub) { return (GenSandiaSubD<Scalar> *) this->subDomain[iSub]; }
  void subConstraintProduct(int iSub, int num_vect, const double* R[], Scalar** V, int trans);
  void findNumGlobNodes();
  void setMpcRhs(int iSub, GenDistrVector<Scalar> &cu);
  void updateMpcRhs(int iSub, GenDistrVector<Scalar> &cu);
  void zeroMpcRhs(int iSub);
  bool setNodeToNodeInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                                 enum face_type *faceEl_type, int *faceEl_ptr, int *faceEl_conn,
                                 double (*faceEl_normal)[3], double (*gap_vector)[3]);
  void setMortarInteractions(int itype, int numNeighb, int *neighb, int *neighbPtr,
                             enum face_type *faceEl_type, int *faceEl_ptr, int *faceEl_conn,
                             int *glbl_surface_id, double *normal_tol, double *tangent_tol);
  void augmentMpclocal();
};

#ifdef _TEMPLATE_FIX_
  #include <Interface.d/SandiaDomain.C>
  #include <Interface.d/SandiaSubD.C>
#endif

#endif

