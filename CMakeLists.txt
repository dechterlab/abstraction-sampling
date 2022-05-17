# CMake file for CVO

cmake_minimum_required(VERSION 2.8)
project(AbstractionSampling)

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
  ./AND-OR
  ./BE
  ./CVO
  ./Problem
  ./Utils
  ./miniSAT
)


# Main executable
add_executable(AbstractionSampling
  ./AbsSampling.cpp
  ./AND-OR/AbsSamplingDFS.cpp
  ./AND-OR/AbsSamplingDFSnonProper.cpp
  ./AND-OR/AbsSamplingWorkspace.cpp
  ./AND-OR/SearchSpaceNodes.cpp
  ./AND-OR/DFSANDcontainer.cpp
  ./BE/MiniBucket.cpp
  ./BE/Bucket.cpp
  ./BE/MBEworkspace.cpp
  ./CVO/Graph.cpp
  ./CVO/Graph_MinFillOrderComputation.cpp
  ./CVO/VariableOrderComputation.cpp
  ./CVO/Graph_RemoveRedundantFillEdges.cpp
  ./Problem/Problem.cpp
  ./Problem/Globals.cpp
  ./Problem/Workspace.cpp
  ./Problem/Function.cpp
  ./Utils/Array.cpp
  ./Utils/AVLtree.cpp
  ./Utils/AVLtreeInt64.cpp
  ./Utils/AVLtreeObjects.cpp
  ./Utils/AVLtreeSimple.cpp
  ./Utils/Heap.cpp
  ./Utils/HeapBasic.cpp
  ./Utils/HeapBasicReductionWithPermutationMapping.cpp
  ./Utils/HierarchicalCluster.cxx
  ./Utils/Mutex.cpp
  ./Utils/MiscUtils.cpp
  ./Utils/FnExecutionThread.cpp
  ./Utils/Sort.cxx
  ./miniSAT/Solver.cpp

)
if (LINK_STATIC)
  SET_TARGET_PROPERTIES(AbstractionSampling PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(AbstractionSampling PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
target_link_libraries(AbstractionSampling ${CMAKE_THREAD_LIBS_INIT})
