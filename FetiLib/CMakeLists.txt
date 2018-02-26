cmake_minimum_required(VERSION 3.5  FATAL_ERROR)

ADD_LIBRARY(FDPSolver "" ConcreteSub.cpp ConcreteSub.h)

TARGET_SOURCES(FDPSolver  PRIVATE FDPSolver.cpp Subdomain.cpp PUBLIC FDPSolver.h DOFInfo.h Subdomain.h)
TARGET_LINK_LIBRARIES(FDPSolver Eigen3::Eigen)