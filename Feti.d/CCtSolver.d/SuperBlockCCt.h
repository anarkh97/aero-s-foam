#ifndef _NEWBLOCKCCT_H_
#define _NEWBLOCKCCT_H_

#include <Feti.d/CCtSolver.d/CCtSolver.h>

template<class Scalar> 
class SuperBlockCCtSolver : public CCtSolver<Scalar>
{
  public:
    SuperBlockCCtSolver(Connectivity *blockToMpc, Connectivity *mpcToMpc, Connectivity *mpcToSub, 
                      Connectivity *mpcToCpu, int numSubsWithMpcs, GenSubDomain<Scalar> **subsWithMpcs,
                      FetiInfo *finfo, FSCommunicator *fetiCom, bool super_flag = true, bool sub_flag = false);
    ~SuperBlockCCtSolver();
    void reSolve(GenDistrVector<Scalar> &v);
    void zeroAll();
    void assemble();
    void factor();

  private:
    FetiInfo *finfo;
    int nMpcBlocks;         // total nb of blocks
    int nMpcBlocksOnMyCPU;
    GenSolver<Scalar> **blockCCtsolver;// array[nMpcBlocks] of pointer on the BlockCCtsolver of myCPU
    Connectivity *mpcToMpc;
    Connectivity *blockToMpc;
    Connectivity *blockToSub;
    Connectivity *mpcToBlock;
    Connectivity **blockMpcToMpc;
    Connectivity *blockToMpcCpu;
    Connectivity *blockToCpu;
    Connectivity *cpuToBlock;
    GenVector<Scalar> **mpcv;
    FSCommPattern<Scalar> *blockCCtPat;
    FSCommPattern<Scalar> *mpcvPat1, *mpcvPat2;
    int myCPU, numCPUs;

    int *nBigBlocksperMPI;                // nb of "big" blocks on myCPU
    int *nSmallBlocksperMPI;              // nb of "small" blocks on myCPU
    Connectivity *MPITosuperBlock;        // for each CPU give the (super) blocks Id they will store & solve
    Connectivity *mpcCpuToBlock;          // for each CPU, give the lmpc block whose have 
                                          //    lmpc CCt contribution from this CPU   
    SimpleNumberer **blockMpcEqNums;      // array[nMpcBlock] of pointer to the EqNum of each blocks created on myCPU 
    int **GlMpcToLlBlkMpcNumMap;          // array[nMpcBlock] of pointer to the lmpc global Id to local block Id mapping 
                                          // of each blocks created on myCPU 
    int nLocAssBlocksOnMyCPU;             // nb of blocks assembled (partially of fully) on myCPU 
    int nTmpAssBlocksOnMyCPU;             // nb of blocks PARTIALLY assembled on myCPU 
    int nExtMpcBlocksOnMyCPU;             // nb of blocks created (i.e. associated with a CCt solver) on myCPU
    ResizeArray<int> *myCPUToLocAssBlocks;// list of blocks (global Id) assembled (partially of fully) on myCPU 
    ResizeArray<int> *myCPUToTmpAssBlocks;// list of blocks (global Id) PARTIALLY assembled on myCPU 
    int *myCPUExtBlockIdArray;            // list of blocks created on myCPU
    Connectivity *LlIdLocAssBlkToMySubs;  // map a locally assembled block (local) Id to the subs of myCPU that contribute to it 
    int *LlIdLocAssBlkToSendId;           // map a locally assembled block (local) Id to their sending (communication) Id
    Connectivity *myCPUBlkToCpuSharedMpcs;// gives the "Cpu-shared" lmpc eqs for each blocks stored & solved on myCPU  
 
    void initialize();
    void makeTmpBlockMpcToMpcConnectivity();
    void makeSubBlocks();
    double estimateBlockCost(compStruct &renumber, double *blockCost, double *blockBandWidth);
    double estimateBlockCost(double *blockCost, double *blockBandWidth);
    void makeSuperBlocks();
    void createSuperBlockCCt(Connectivity *mpcToSub);
    void makeBlockCCtCommPattern();
    void setBlockCCtCommSize(int iBlock);
    void factorSmallBlockCCtsolver(int iBlock, int *blockOffset);
    void factorSparseBigBlockCCtsolver(int iBlock);
    void createBlockMpcToMpcConnectivity();
    void createOneBlockMpcToMpcConnectivity(int IBlock);
    void makeGlobalMpcToLocalBlkMpcNumMap();
    void makeOneGlobalMpcToLocalBlkMpcNumMap(int IBlock);
    void createBlockCCtsolver();
    void createOneBlockCCtsolver(int IBlock);
    void assembleOneBlockCCtsolver(int IBlock);
    void zeroOneBlockCCtsolver(int IBlock);
    void deleteGlobalMpcToLocalBlkMpcNumMap();
    void deleteOneGlobalMpcToLocalBlkMpcNumMap(int IBlock);
    void deleteBlockMpcToMpcConnectivity();
    void deleteMyCPUTmpBlockCCtsolver();
    void deleteOneBlockMpcToMpcConnectivity(int IBlock);
    void deleteBlockMpcEqNums();
    void deleteOneBlockMpcEqNums(int IBlock, int *gBlockIdArray);
    void deleteBlockCCtsolver();
    void deleteOneBlockCCtsolver(int IBlock, int *gBlockIdArray);
    void deleteMpcVector();
    void initBlockMpcResidual();
    void solveBlockCCt(GenDistrVector<Scalar> &v);
    void extractBlockMpcResidual(GenDistrVector<Scalar> &v);
    void insertBlockMpcResidual(int iSub, GenDistrVector<Scalar> &v);
    void initOneBlockMpcResidual(int IBlock);
    void solveOneBlockCCt(int IBlock, GenDistrVector<Scalar> &v);
    void extractOneBlockMpcResidual(int IBlock, GenDistrVector<Scalar> &v);
    void sendOneBlockCCtsolver(int IBlock);
    void recOneBlockCCtsolver(int IBlock);
    void sendBlockMpcResidualBeforeSolve();
    void recBlockMpcResidualBeforeSolve();
    void sendBlockMpcResidualAfterSolve();
    void recBlockMpcResidualAfterSolve();
    void sendOneBlockMpcResidualBeforeSolve(int IBlock);
    void recOneBlockMpcResidualBeforeSolve(int IBlock);
    void sendOneBlockMpcResidualAfterSolve(int IBlock);
    void recOneBlockMpcResidualAfterSolve(int IBlock);
    void makeBlkToCpuSharedMpcsMap();
};

#ifdef _TEMPLATE_FIX_
  #include<Feti.d/CCtSolver.d/SuperBlockCCt.C>
#endif

#endif
