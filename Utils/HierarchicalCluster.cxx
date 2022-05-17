/*
	This module contains hierchical cluseter function's declarations.

	Bobak Pezeshki, April 2022.
*/

#include "HierarchicalCluster.hxx"
#include "LowerTriangularMatrix.hxx"
#include "Problem.hxx"
#include <cmath>
#include <utility>
#include <queue>
#include <set>
#include <iterator>
#include <cassert>


struct GreaterCmpWithNextCluster{
    bool operator()(std::pair<uint32_t,double> a, std::pair<uint32_t,double> b) const{
        return a.second > b.second;
    }
};

void SortedLinearMinVarianceHierchicalClustering(std::vector<uint32_t>& exclusiveEndPoints, void *objarray[], uint32_t len, uint32_t numClusters, ValueToClusterByFn valFn){
    if(len <= numClusters){
        exclusiveEndPoints.reserve(len);
        for(uint32_t idx = 1; idx <= len; ++idx){
            exclusiveEndPoints.push_back(idx);
        }
    }
    else{
        //initialization
        exclusiveEndPoints.reserve(numClusters);
        std::set<uint32_t> sorted_clusters;
        LowerTriangularMatrix<double> distMatrix(len);
        std::vector<double> adjacentDist(len);
        std::vector<uint32_t> clusterSizes(len,1);
        std::priority_queue<std::pair<uint32_t,double>, std::vector<std::pair<uint32_t,double>>, GreaterCmpWithNextCluster > minDist;
        uint32_t clusterAssignment[len];
        for(int32_t idx = len - 2; idx >= 0; --idx){
            sorted_clusters.emplace(idx);
            clusterAssignment[idx] = idx;
            for(uint32_t idxx = idx+1; idxx < len; ++idxx){
                double logLinearDist;
                // LOG_OF_SUB_OF_TWO_NUMBERS_GIVEN_AS_LOGS(logLinearDist, valFn(objarray[idx]), valFn(objarray[idxx]))
                LOG_OF_SUB_OF_TWO_NUMBERS_GIVEN_AS_LOGS(logLinearDist, (*valFn)(objarray[idx]), (*valFn)(objarray[idxx]))
                distMatrix(idx,idxx) = 2 * logLinearDist;
            }
            adjacentDist[idx] = distMatrix(idx,idx+1);
            minDist.emplace(idx, adjacentDist[idx]);
        }
        //last cluster will have no distMatrix or adjacencyDist distance to a successive cluster (and never be merged into)
        uint32_t idx = len -1;
        sorted_clusters.emplace(idx);
        clusterAssignment[idx] = idx;
        
        //merge clusters
        uint32_t numRounds = len - numClusters;
        for(uint32_t round = 0; round < numRounds; ++round){
            uint32_t idx1;
            double dist1;
            do{
                idx1 = minDist.top().first;
                dist1 = minDist.top().second;
                minDist.pop();
            } while(dist1 != adjacentDist[idx1]);
            std::set<uint32_t>::iterator idx1_fwd_it = sorted_clusters.find(idx1);
            std::reverse_iterator<std::set<uint32_t>::iterator> idx0_rev_it(idx1_fwd_it);
            std::set<uint32_t>::iterator idx2_fwd_it = std::next(idx1_fwd_it);  //idx2 = cluster immediately after idx1
            std::set<uint32_t>::iterator idx3_fwd_it = std::next(idx2_fwd_it);  //idx3 = cluster immediately after idx2
            uint32_t idx2 = *( idx2_fwd_it );
            //uint32_t idx2end = clusterExclusiveEndpoint[idx2]; //=len if there is no other cluster after cluster idx2
            //clusterExclusiveEndpoint[idx1] = idx2end;
            sorted_clusters.erase(idx2_fwd_it);
            adjacentDist[idx2] = NAN; //implicitly markes cluster idx2 as removed
            for(std::reverse_iterator<std::set<uint32_t>::iterator> rev_it = idx0_rev_it; rev_it != sorted_clusters.rend(); ++rev_it){
                //update distances from all remaining sorted_clusters < idx1 to idx1
                uint32_t idx = *rev_it;
                distMatrix(idx,idx1) = linearWardsDistanceUpdate(
                    clusterSizes[idx1],//n_i, 
                    clusterSizes[idx2],//n_j, 
                    clusterSizes[idx],//n_k
                    distMatrix.val(idx,idx1),//D_ki, 
                    distMatrix.val(idx,idx2),//D_kj, 
                    distMatrix.val(idx1,idx2)//D_ij
                );
            }
            for(std::set<uint32_t>::iterator fwd_it = idx3_fwd_it; fwd_it != sorted_clusters.end(); ++fwd_it){
                //update distances from cluster idx1 to all remaining sorted_clusters > idx1
                uint32_t idx = *fwd_it;
                distMatrix(idx,idx1) = linearWardsDistanceUpdate(
                    clusterSizes[idx1],//n_i, 
                    clusterSizes[idx2],//n_j, 
                    clusterSizes[idx],//n_k
                    distMatrix.val(idx1,idx),//D_ik, 
                    distMatrix.val(idx2,idx),//D_jk, 
                    distMatrix.val(idx1,idx2)//D_ij
                );
            }
            clusterSizes[idx1]+=clusterSizes[idx2];
            //update adjacentDist[idx1]
            if (idx3_fwd_it != sorted_clusters.end()){
                uint32_t idx3 = *idx3_fwd_it;
                adjacentDist[idx1] = distMatrix(idx1,idx3);
                minDist.emplace(idx1, adjacentDist[idx1]);
            }
            else{
                adjacentDist[idx1] = NAN; //ensures will never be valid choice from the priority queue
            }
            //update adjacentDist[idx0]
            if (idx0_rev_it != sorted_clusters.rend()){
                uint32_t idx0 = *idx0_rev_it;
                adjacentDist[idx0] = distMatrix(idx0,idx1);
                minDist.emplace(idx0, adjacentDist[idx0]);
            }
        }

        std::set<uint32_t>::iterator fwd_it_begin = sorted_clusters.begin();
        std::set<uint32_t>::iterator fwd_it_end = ++fwd_it_begin;
        for (fwd_it_begin, fwd_it_end; fwd_it_end != sorted_clusters.end(); fwd_it_begin = fwd_it_end, ++fwd_it_end){
            exclusiveEndPoints.emplace_back(*fwd_it_end);
        }
        exclusiveEndPoints.emplace_back(len);
        assert(sorted_clusters.size() == numClusters);
    }
}

double linearLanceWilliamsDistanceUpdate(double a_i, double D_ik, double a_j, double D_jk, double B, double D_ij, double g){
    return a_i * D_ik + a_j * D_jk + B * D_ij - g * std::abs(D_ik - D_jk);
}

double linearWardsDistanceUpdate(uint32_t n_i, uint32_t n_j, uint32_t n_k, double D_ik, double D_jk, double D_ij){
    double n_ijk = n_i + n_j + n_k;
    return (n_i + n_k)/n_ijk * D_ik + (n_j + n_k)/n_ijk * D_jk - n_k/n_ijk * D_ij;
}
