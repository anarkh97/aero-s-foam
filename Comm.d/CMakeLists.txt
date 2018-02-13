cmake_minimum_required(VERSION 3.5  FATAL_ERROR)
ADD_LIBRARY(Comm BaseCommunicator.cpp BaseCommunicator.h OpaqueHandle.cpp OpaqueHandle.h MPICompatTraits.h)
#
#SET(INC_DIR ..)
#
#INCLUDE_DIRECTORIES( ${INC_DIR} . )

TARGET_SOURCES(Comm
        PRIVATE CommunicatorCore.C
        Communicator.h
        )

ADD_EXECUTABLE(comm_test CommTest.cpp)
#TARGET_SOURCES(comm_test PUBLIC )

TARGET_LINK_LIBRARIES(comm_test Comm ${MPI_LIBRARIES})