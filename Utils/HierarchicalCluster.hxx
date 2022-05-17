/*
	This module contains hierchical cluseter function's declarations.

	Bobak Pezeshki, April 2022.
*/

#ifndef HIERARCHICAL_CLUSTER_HXX_INCLUDED
#define HIERARCHICAL_CLUSTER_HXX_INCLUDED

#include <inttypes.h>
#include <vector>

#ifdef LINUX
typedef int64_t INT64 ;
#else 
typedef __int64 INT64 ;
#endif

// typedef std::function<double(void *)> ValueToClusterByFn;
// void SortedLinearMinVarianceHierchicalClustering(std::vector<uint32_t>& exclusiveEndPoints, void *objarray[], uint32_t len, uint32_t numClusters, ValueToClusterByFn valFn);


typedef double (*ValueToClusterByFn)(void *Obj1) ;
void SortedLinearMinVarianceHierchicalClustering(std::vector<uint32_t>& exclusiveEndPoints, void *objarray[], uint32_t len, uint32_t numClusters, ValueToClusterByFn valFn);
double linearLanceWilliamsDistanceUpdate(double a_i, double D_ik, double a_j, double D_jk, double B, double D_ij, double g);
double linearWardsDistanceUpdate(uint32_t n_i, uint32_t n_j, uint32_t n_k, double D_ik, double D_jk, double D_ij);

#endif // SORT_HXX_INCLUDED
