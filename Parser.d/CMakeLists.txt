# search bison
#MACRO(FIND_BISON)
#    IF(NOT BISON_EXECUTABLE)
#        FIND_PROGRAM(BISON_EXECUTABLE bison)
#        IF (NOT BISON_EXECUTABLE)
#          MESSAGE(FATAL_ERROR "flex not found - aborting")
#        ENDIF (NOT BISON_EXECUTABLE)
#    ENDIF(NOT BISON_EXECUTABLE)
#ENDMACRO(FIND_BISON)
#
#MACRO(ADD_BISON_FILES _sources )
#    FIND_BISON()
#
#    FOREACH (_current_FILE ${ARGN})
#      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
#      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)
#
#      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.tab.C)
#
#      ADD_CUSTOM_COMMAND(
#         OUTPUT ${_out}
#         COMMAND ${BISON_EXECUTABLE}
#         ARGS
#         -d
#         -v
#         -p yy
#         -o${_out}
#         ${_in}
#         DEPENDS ${_in}
#      )
#
#      SET(${_sources} ${${_sources}} ${_out} )
#   ENDFOREACH (_current_FILE)
#ENDMACRO(ADD_BISON_FILES)


# flex a .ll file

# search flex
#MACRO(FIND_FLEX)
#    IF(NOT FLEX_EXECUTABLE)
#        FIND_PROGRAM(FLEX_EXECUTABLE flex)
#        IF (NOT FLEX_EXECUTABLE)
#          MESSAGE(FATAL_ERROR "flex not found - aborting")
#        ENDIF (NOT FLEX_EXECUTABLE)
#    ENDIF(NOT FLEX_EXECUTABLE)
#ENDMACRO(FIND_FLEX)
#
#MACRO(ADD_FLEX_FILES _sources )
#    FIND_FLEX()
#
#    FOREACH (_current_FILE ${ARGN})
#      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
#      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)
#
#      SET(_out ${CMAKE_CURRENT_BINARY_DIR}/flex_${_basename}.cpp)
#
#      ADD_CUSTOM_COMMAND(
#         OUTPUT ${_out}
#         COMMAND ${FLEX_EXECUTABLE}
#         ARGS
#         -Pyy
#         -o${_out}
#         ${_in}
#         DEPENDS ${_in}
#      )
#
#      SET(${_sources} ${${_sources}} ${_out} )
#   ENDFOREACH (_current_FILE)
#ENDMACRO(ADD_FLEX_FILES)
#
#ADD_FLEX_FILES(MyLib lexer.l)
#ADD_BISON_FILES(MyLib p.y)
#
#ADD_LIBRARY(Parser AuxDefs.C ${MyLib})
#SET_TARGET_PROPERTIES(Parser
#                PROPERTIES COMPILE_FLAGS -fPIC)

find_package(BISON)
find_package(FLEX)
BISON_TARGET(MyParser p.y ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp)
FLEX_TARGET(MyScanner lexer.l  ${CMAKE_CURRENT_BINARY_DIR}/lexer.cpp)
ADD_FLEX_BISON_DEPENDENCY(MyScanner MyParser)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
ADD_LIBRARY(Parser AuxDefs.C ${BISON_MyParser_OUTPUTS} ${FLEX_MyScanner_OUTPUTS})
SET_TARGET_PROPERTIES(Parser PROPERTIES COMPILE_FLAGS -fPIC) 
