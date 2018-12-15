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
        MultiDomainSolver.h MultiDomainSolver.C
		)
if(USE_CHOLMOD)
	message("Looking for CHOLMOD")
	find_package(Cholmod)
	message("Cholmod include: ${CHOLMOD_INCLUDE_DIR}")
	message("Cholmod libraries: ${CHOLMOD_LIBRARIES}")
	if (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
		message("Will use CholMod")
		target_sources(solver
				PRIVATE
				CholmodSolver.cpp
				CholmodImp.h CholmodImp.cpp
				PUBLIC
				CholmodSolver.h)
		target_compile_definitions(solver PUBLIC WITH_CHOLMOD)
		target_include_directories(solver PRIVATE ${CHOLMOD_INCLUDE_DIR})
	endif (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
endif(USE_CHOLMOD)

TARGET_LINK_LIBRARIES(solver Eigen3::Eigen)

#SET(INC_DIR ..)

#INCLUDE_DIRECTORIES( ${INC_DIR} )

