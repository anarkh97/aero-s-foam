ADD_LIBRARY(Feti
        CoarseSetCore.C
        CGOrthoSet.C CGOrthoSet.h
        DistrVectorCore.C
        FetiCore.C
        FetiOpCore.C
        NLFetiCore.C
        CCtSolver.d/BlockCCt.C
        FetiSub.h)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
