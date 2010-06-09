/* DEPRECATED, we now use SingleDomainDynam in DynamDescr.h instead
#ifndef _TEMP_DESCR_H_
#define _TEMP_DESCR_H_

#include <Driver.d/Domain.h>
#include <Driver.d/TempProbType.h>
#include <Driver.d/Dynam.h>
#include <Hetero.d/FlExchange.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class DynOps;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenFSFullMatrix;
typedef GenFSFullMatrix<double> FSFullMatrix;


class SDTempDynamPostProcessor {
    Domain *domain;
    double *bcx;
  public:
    SDTempDynamPostProcessor(Domain *d, double *_bcx) 
	{ domain = d; bcx = _bcx; }
    void tempdynamOutput(int,DynamMat&,Vector&,TempState<Vector>&);
};

class SingleDomainTemp {
    Domain *domain;
    int    *bc;
    double *bcx;
    int    *nodeOrder;
    SparseMatrix *kuc;

    //For temproject
    FSFullMatrix  *X;    // pre-calculated projector
    double *Rmem;        // global Zero Energy modes 
    int numR;            // number of Zero Energy Modes

    int maxmode;      // number of eigenvectors for modal decomposition analysis
    double **eigmodes;
    double **tPhiM;   // Stores Transpose(Phi)*M

 public:
    SingleDomainTemp(Domain *d) { domain = d; }
    int* boundary() { return bc;}
    double* boundaryValue() { return bcx;}

    DynamMat buildOps(double, double, double);
//    DynamMat buildDynamicOps(double, double, double);
    Solver *getSolver();

    SDTempDynamPostProcessor *getPostProcessor();
    int solVecInfo();
    int dbcVecInfo();
    int getTimeIntegration();
    int getAeroheatFlag();
    int getThermohFlag();
    int getHzemFlag();
    int getZEMFlag();

    void temptrProject(Vector &f);
    void tempProject(Vector &v);
    void tempprojector_prep(Rbm *R, SparseMatrix *M);
    void getTempTimes(double &dtemp, double &tmax);
    void getFluiddTime(double &dtfluid, double &tmaxfluid);
    void getTempNeum(double &epsiln);
    void tempInitState(TempState<Vector> &);
    void computeExtForce(Vector &, double t, int tIndex, Vector &);
    void preProcess();
    void aeroHeatPreProcess(Vector& d_n, Vector& v_n, Vector& v_p);
    void thermohPreProcess(Vector& d_n, Vector& v_n, Vector& v_p);
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol);
    void getQuasiStaticParameters(double &maxVel, double &qsbeta);
    void getInternalForce(Vector&, Vector&);
    void getInitialTime(int &tIndex, double &initialTime);

    int getModeDecompFlag();
    void modeDecompPreProcess(SparseMatrix *M);
    void modeDecomp(double t, int tIndex, Vector& d_n);

    FlExchanger *flExchanger;
    int cmdComHeat(int cmdFlag) { return flExchanger->cmdComHeat(cmdFlag); }

};


#endif
*/
