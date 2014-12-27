# -*-  mode: CMAKE; cmake-tab-width: 4 -*-
#

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
ADD_LIBRARY(glps parser scanner elemtree.c glps.c glpserror.c glpsutil.cc)
#SET(CMAKE_CXX_FLAGS "-ggdb3")
#SET(CMAKE_C_FLAGS "-ggdb3")

# parse the full version number from zlib.h and include in ZLIB_FULL_VERSION
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/glps.h _glps_h_contents)
string(REGEX REPLACE ".*#define[ \t]+GLPS_VERSION[ \t]+\"([0-9A-Za-z.]+)\".*"
    "\\1" GLPS_FULL_VERSION ${_glps_h_contents})

#MESSAGE(STATUS "------------")
#MESSAGE(STATUS "glps version: ${GLPS_FULL_VERSION}")
#CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/glps/glps_config.h.in"
#    "${PROJECT_BINARY_DIR}/glps/glps_config.h"
#)

## Create target for the parser

FIND_PROGRAM(LEX_EXE flex)
IF(LEX_EXE STREQUAL "LEX_EXE-NOTFOUND")
    MESSAGE(FATAL_ERROR "ERROR: flex not found!")
ENDIF(LEX_EXE STREQUAL "LEX_EXE-NOTFOUND")

FIND_PROGRAM(YACC_EXE bison )
IF(YACC_EXE STREQUAL "YACC_EXE-NOTFOUND")
  MESSAGE(FATAL_ERROR "ERROR: yacc not found!")
ENDIF(YACC_EXE STREQUAL "YACC_EXE-NOTFOUND")

MESSAGE(STATUS "Found lex : ${LEX_EXE}")
MESSAGE(STATUS "Found yacc: ${YACC_EXE}")

# Create custom command for flex/lex (note the outputs)
ADD_CUSTOM_COMMAND(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/scanner.l
  COMMAND ${LEX_EXE} 
  ARGS -i -o${CMAKE_CURRENT_BINARY_DIR}/scanner.c ${CMAKE_CURRENT_SOURCE_DIR}/scanner.l
  TARGET dak
  OUTPUTS ${CMAKE_CURRENT_BINARY_DIR}/scanner.c)

# Create custom command for bison/yacc (note the DEPENDS)
ADD_CUSTOM_COMMAND(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/parser.y
  COMMAND ${YACC_EXE}
  ARGS -t -v -d -y ${CMAKE_CURRENT_SOURCE_DIR}/parser.y -o ${CMAKE_CURRENT_BINARY_DIR}/parser.c
  TARGET dak
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/scanner.c
  OUTPUTS ${CMAKE_CURRENT_BINARY_DIR}/parser.c)
#endif(WIN32)
# Add parser.c to the list of sources
#SET(dak_SRCS ${dak_SRCS} ${dak_BINARY_DIR}/src/parser.c beamline.cc latio.cc)

# Since parser.c does not exists yet when cmake is run, mark
# it as generated
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/parser.c GENERATED)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/parser.h GENERATED)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/scanner.c GENERATED)

#INSTALL(TARGETS glps DESTINATION lib)
#INSTALL(FILES glps.h glpserror.h glpsutil.h
#    DESTINATION include/glps)

#IF(MSVC)
#    SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/src/parser.c 
#	COMPILE_FLAGS "/Tp")
#ENDIF(MSVC)

## Include binary directory to include lexer.c in parser.c
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/src)

