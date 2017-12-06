ADD_LIBRARY(Feti
        CoarseSetCore.C
        CoarseSet.C
        CoarseSet.h
        CGOrthoSet.C CGOrthoSet.h
        DistrVectorCore.C
        FetiCore.C
        FetiOpCore.C
        FetiDP.C
        Feti.h
        NLFetiCore.C
        CCtSolver.d/BlockCCt.C
        FetiSub.h FETISolver.cpp FETISolver.h FetiSub.cpp)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
