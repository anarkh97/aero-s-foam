ADD_LIBRARY(solver
        Solver.h
        DSCsolver.C
        MumpsCore.C
        Rbm.C
        MappedAssembledSolver.C
        SolverCore.C
        SpoolesCore.C
        GoldfarbIdnani.C
        SolverFactoryCore.C
        MultiDomainSolver.h MultiDomainSolver.C)

#SET(INC_DIR ..)

#INCLUDE_DIRECTORIES( ${INC_DIR} )

