# CMake file for CVO

cmake_minimum_required(VERSION 2.8)
project(ARP-v3.0.6)

# to enable static linking
option(LINK_STATIC "Link binary statically" ON)

#add_definitions(-DDEFINE_PACE16_MAIN_FN)

if(WIN32)
  add_definitions(-DWINDOWS)
else()
  add_definitions(-DLINUX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -w -fpermissive -fPIC")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -fPIC")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w -fpermissive -fPIC")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -fPIC")
  set(CMAKE_CXX_FLAGS "-ggdb")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG /Od /J")
  set(CMAKE_C_FLAGS_DEBUG "-DDEBUG /Od /J")
  set(CMAKE_CXX_FLAGS_RELEASE "/Ox /J")
  set(CMAKE_C_FLAGS_RELEASE "/Ox /J")
endif()


if(NOT WIN32)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
  else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif()

if(LINK_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

find_package(Threads)

include_directories(
  ./
  ../../../ARP/ProgramVersions/ARP-v3.0.6
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR
  ../../../ARP/ProgramVersions/ARP-v3.0.6/BE
  ../../../ARP/ProgramVersions/ARP-v3.0.6/CVO
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Problem
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils
  ../../../miniSAT
)


# Main executable
add_executable(ARP-v3.0.6
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AbsSampling.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR/AbsSamplingDFS.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR/AbsSamplingDFSnonProper.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR/AbsSamplingWorkspace.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR/SearchSpaceNodes.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/AND-OR/DFSANDcontainer.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/BE/MiniBucket.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/BE/Bucket.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/BE/MBEworkspace.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/CVO/Graph.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/CVO/Graph_MinFillOrderComputation.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/CVO/VariableOrderComputation.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/CVO/Graph_RemoveRedundantFillEdges.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Problem/Problem.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Problem/Globals.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Problem/Workspace.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Problem/Function.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/Array.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/AVLtree.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/AVLtreeInt64.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/AVLtreeObjects.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/AVLtreeSimple.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/Heap.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/HeapBasic.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/HeapBasicReductionWithPermutationMapping.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/Mutex.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/MiscUtils.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/FnExecutionThread.cpp
  ../../../ARP/ProgramVersions/ARP-v3.0.6/Utils/Sort.cxx
  ../../../miniSAT/Solver.cpp
)
if (LINK_STATIC)
  SET_TARGET_PROPERTIES(ARP-v3.0.6 PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(ARP-v3.0.6 PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
target_link_libraries(ARP-v3.0.6 ${CMAKE_THREAD_LIBS_INIT})
