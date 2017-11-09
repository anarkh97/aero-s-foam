#ifndef _BLOCKCCT_H_
#define _BLOCKCCT_H_
#include <vector>
#include <Feti.d/CCtSolver.d/CCtSolver.h>

template<class Scalar>
class BlockCCtSolver : public CCtSolver<Scalar>
{
  public:
    BlockCCtSolver(Connectivity *blockToMpc, Connectivity *mpcToMpc, Connectivity *mpcToSub,
                   Connectivity *mpcToCpu, int numSubsWithMpcs, GenSubDomain<Scalar> **subsWithMpcs,
                   int *subMap, FetiInfo *finfo, FSCommunicator *fetiCom);
    ~BlockCCtSolver();
    void reSolve(GenDistrVector<Scalar> &v);
    void zeroAll();
    void assemble();
    void factor();

  private:
    FetiInfo *finfo;
    int nMpcBlocks;         // total nb of blocks
    int nMpcBlocksOnMyCPU;  // nb of blocks stored & solved on myCPU
    int *subMap;            // map from global subdomain index to subsWithMpcs index
    std::vector<GenSolver<Scalar> *> blockCCtsolver;// array[nMpcBlocks] of pointer on the BlockCCtsolver of myCPU
    std::vector<GenSparseMatrix<Scalar> *> blockCCtsparse;
    std::vector<SimpleNumberer *> blockMpcEqNums;
    Connectivity *blockToMpc;
    Connectivity *blockToSub;
    Connectivity *mpcToBlock;
    std::vector<Connectivity *> blockMpcToMpc;
    Connectivity *blockToMpcCpu;   //HB: for each block, give the CPUs that have a 
                                   //    lmpc CCt contribution to it
    Connectivity *blockToCpu;
    Connectivity *cpuToBlock;
    mutable std::vector<GenVector<Scalar> *> mpcv; // This make this solver non re-entrant.
    FSCommPattern<Scalar> *blockCCtPat;
    FSCommPattern<Scalar> *mpcvPat1, *mpcvPat2;
    int myCPU, numCPUs;

    void createBlockMpcToMpcConnectivity(int iBlock, Connectivity *mpcToMpc);
    void deleteBlockMpcToMpcConnectivity(int iBlock);
    void createBlockCCtsolver(int iBlock);
    void deleteBlockCCtsolver(int iBlock);
    void assembleBlockCCtsolver(int i);
    void initBlockMpcResidual(int iBlock) const; // Modifies mutable mpcv
    void makeBlockCCtCommPattern();
    void setBlockCCtCommSize(int iBlock);
    void sendBlockCCtsolver(int iBlock);
    void recBlockCCtsolver(int iBlock);
    void sendBlockMpcResidualBeforeSolve(int iBlock);
    void recBlockMpcResidualBeforeSolve(int iBlock);
    void sendBlockMpcResidualAfterSolve(int iBlock);
    void recBlockMpcResidualAfterSolve(int iBlock);
    void distributeMpcBlocks();
    void factorBlockCCtsolver(int iBlock);
    void extractSolveAndInsertBlockMpcResidual(GenDistrVector<Scalar> &v);
    void solveBlockCCt(int iBlock, GenDistrVector<Scalar> &v) const;
    void extractBlockMpcResidual(int iBlock, GenDistrVector<Scalar> &v) const;
    void insertBlockMpcResidual(int i, GenDistrVector<Scalar> &v);
    void zeroBlockCCtsolver(int i);
};

#endif
