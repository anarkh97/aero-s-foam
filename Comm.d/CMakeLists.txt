cmake_minimum_required(VERSION 3.5  FATAL_ERROR)
ADD_LIBRARY(Comm OpaqueCommunicator.cpp OpaqueCommunicator.h OpaqueHandle.cpp OpaqueHandle.h)
#
#SET(INC_DIR ..)
#
#INCLUDE_DIRECTORIES( ${INC_DIR} . )

TARGET_SOURCES(Comm
        PRIVATE CommunicatorCore.C
        Communicator.h
        )

