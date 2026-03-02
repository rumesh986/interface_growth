# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_dim subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(one_dim C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

# Add a block of this type for your own executable
# oomph_add_executable(
#   NAME erf2d18 # <--- name of the executable
#   SOURCES erf2d18.cc includes.h two_phase_free_boundary_mesh.h # <--- source code(s)
#   LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic) # <-- oomph-lib
#                                                          # libraries required
#                                                          # (order matters:
#                                                          # specific ones go
#                                                          # first; more general
#                                                          # ones at the end)

function(cust_add_executable)
  set(PREFIX ARG)
  set(FLAGS JUST_IN_CASE)
  set(SINGLE_VALUE_ARGS TARGET X_ORDER T_ORDER)
  set(MULTI_VALUE_ARGS SOURCES)

  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  oomph_add_executable(
    NAME ${ARG_TARGET}_${ARG_X_ORDER}_${ARG_T_ORDER} # <--- name of the executable
    SOURCES ${ARG_SOURCES} # <--- source code(s)
    LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic
    CXX_OPTIONS -Wall -Werror
    CXX_DEFINITIONS RUN_SCRIPT X_ORDER=${ARG_X_ORDER} T_ORDER=${ARG_T_ORDER}
  )
endfunction(cust_add_executable)

function(cust_add_executables)
  set(PREFIX ARG)
  set(FLAGS JUST_IN_CASE)
  set(SINGLE_VALUE_ARGS TARGET)
  set(MULTI_VALUE_ARGS SOURCES)

  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  cust_add_executable(
    TARGET ${ARG_TARGET}
    SOURCES ${ARG_SOURCES}
    X_ORDER 2
    T_ORDER 1
  )

  cust_add_executable(
    TARGET ${ARG_TARGET}
    SOURCES ${ARG_SOURCES}
    X_ORDER 3
    T_ORDER 1
  )

  cust_add_executable(
    TARGET ${ARG_TARGET}
    SOURCES ${ARG_SOURCES}
    X_ORDER 2
    T_ORDER 2
  )

  cust_add_executable(
    TARGET ${ARG_TARGET}
    SOURCES ${ARG_SOURCES}
    X_ORDER 3
    T_ORDER 2
  )
endfunction(cust_add_executables)


cust_add_executables(
  TARGET one_dim1
  SOURCES one_dim1.cc two_phase_free_boundary_mesh.h
)

cust_add_executables(
  TARGET two_dim1
  SOURCES two_dim1.cc two_phase_free_boundary_mesh.h
)

message(VERBOSE "Leaving one_dim subdirectory")
# ------------------------------------------------------------------------------
