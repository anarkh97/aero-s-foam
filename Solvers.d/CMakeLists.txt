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
		CholmodSolver.cpp CholmodSolver.h
		)
if(USE_CHOLMOD)
	message("Looking for CHOLMOD")
	find_package(Cholmod)
	message("Cholmod include: ${CHOLMOD_INCLUDE_DIR}")
	message("Cholmod libraries: ${CHOLMOD_LIBRARIES}")
	if (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
		message("Will use CholMod and ${CHOLMOD_FOUND}")
		# pass whether we use choldmod or not to the base CMakeLists.
		set (CHOLMOD_FOUND ${CHOLMOD_FOUND} PARENT_SCOPE)
		target_sources(solver
				PRIVATE
				CholmodImp.h CholmodImp.cpp
				)
		target_compile_definitions(solver PUBLIC WITH_CHOLMOD)
		target_include_directories(solver PRIVATE ${CHOLMOD_INCLUDE_DIR})
	endif (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
endif(USE_CHOLMOD)

TARGET_LINK_LIBRARIES(solver Eigen3::Eigen)

#SET(INC_DIR ..)

#INCLUDE_DIRECTORIES( ${INC_DIR} )

