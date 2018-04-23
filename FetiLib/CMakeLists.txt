cmake_minimum_required(VERSION 3.5  FATAL_ERROR)

ADD_LIBRARY(FDPSolver ConcreteSub.cpp ConcreteSub.h SparseMatrix.cpp SharedNodes.cpp Details/DOFManip.cpp
		Details/SubImpl.cpp Details/SubImpl.h LibTest.cpp)

TARGET_SOURCES(FDPSolver
		PRIVATE FDPSolver.cpp Subdomain.cpp Details/DOFManip.h
		PUBLIC FDPSolver.h DOFInfo.h Subdomain.h
		VectorReference.h SparseMatrix.h Types.h SharedNodes.h
		)
TARGET_LINK_LIBRARIES(FDPSolver Eigen3::Eigen)

ADD_EXECUTABLE(fetiTest LibTest.cpp)
TARGET_LINK_LIBRARIES(fetiTest Comm)