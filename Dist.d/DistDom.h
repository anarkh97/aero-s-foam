#ifndef _DIST_DOM_
#define _DIST_DOM_

#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Utils.d/MyComplex.h>
#include <Driver.d/DecDomain.h>

template<class Scalar>
class GenDistrDomain : virtual public GenDecDomain<Scalar> 
{
  private:
    DistrInfo masterInfo;         // used to create master dist vecs
    int **masterFlag;             // master flag for distributed printout
    int *numFlags;
    int *nodeOffsets;             // offsets for distributed writing
    int *elemNodeOffsets;
    int *elemOffsets;
    int *numRes;		  // number of solutions
    FSCommPattern<Scalar> *nodePat;
    int x;
    DistVec<Scalar> *masterStress;
  public:
    GenDistrDomain(Domain *d);
    virtual ~GenDistrDomain();

    void postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f, double eigV = 0.0,
                        GenDistrVector<Scalar> *aeroF = 0, int x = 0, GenMDDynamMat<Scalar> *dynOps = 0,
                        SysState<GenDistrVector<Scalar> > *distState = 0, int ndflag = 0); 
    void postProcessing(DistrGeomState *u, Corotator ***, double x = 0, SysState<GenDistrVector<Scalar> > *distState = 0,
                        GenDistrVector<Scalar> *aeroF = 0, DistrGeomState *refState = 0);
    virtual void forceContinuity(GenDistrVector<Scalar> &u);

    void setsizeSfemStress(int fileNumber);  // YYY DG implementation incomplete: do the element stresses 
    int getsizeSfemStress() { return this->sizeSfemStress; } // YYY DG for both node- and element-based ?
    Scalar * getSfemStress(int fileNumber); // YYY DG implementation incomplete: do the element stresses 
    void updateSfemStress(Scalar* str, int fileNumber);

  private:
    void initialize();
    void initPostPro();
    void makeNodePat();
    void makeMasterInfo();
    void createMasterFlag();
    void createOutputOffsets();
    void getPrimal(DistSVec<Scalar, 11> &disps, DistSVec<Scalar, 11> &masterDisps,
                   double time, int x, int fileNumber, int ndof, int startdof);//DofSet::max_known_nonL_dof
    void getAeroForceScalar(DistSVec<Scalar, 6> &aerof, DistSVec<Scalar, 6> &masterAeroF,
                            double time, int x, int fileNumber, int dof);
//    void getStressStrain(GenDistrVector<Scalar> &, double, int, int, int);
    void getStressStrain(GenDistrVector<Scalar> &u, double time, int x, int fileNumber, int Findex, int printFlag=0);
    void getStressStrain(GenDistrVector<Scalar> &sol, int fileNumber, int stressIndex, double time, int printFlag=0)  // To keep the same arguments as DecDomain
            { getStressStrain(sol, time, 0, fileNumber, stressIndex, printFlag);} // YYY DG assuming x=0, but what is x ?
    void getElementStressStrain(GenDistrVector<Scalar> &, double, int, int, int, int printFlag=0); // YYY DG Implement printFlag
    void getPrincipalStress(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementPrincipalStress(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementForce(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementAttr(int, int, double);
    void getStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                         int x, int fileNumber, int Findex, DistrGeomState *refState);
    void getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                int iter, int fileNumber, int Findex, DistrGeomState *refState);
    void getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                            int x, int fileNumber, int strIndex, DistrGeomState *refState);
    void getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                                   int x, int fileNumber, int strIndex, DistrGeomState *refState);
    void unify(DistSVec<Scalar, 11> &vec);

};

#ifdef _TEMPLATE_FIX_
  #include <Dist.d/DistDom.C>
#endif

#endif
