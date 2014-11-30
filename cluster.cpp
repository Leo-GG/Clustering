/**
 * @file cluster.cpp
 * @brief Implementation of methods for Cluster class
 *
 * This file implements the calcCentroid method, which assigns the cluster
 * centroid and the cluster radius. The centroid is defined as the element
 * whose sum of distances to all the other members is the smallest.
 */


#include "node.h"
#include "cluster.h"
#include <iostream>
#include <fstream>


void Cluster::calcCentroid(vector< vector<float> > normScores)
{
        float distToNeighbors;
        float minDistanceSum;
        float currentRadius;
        for (int i=0;i<members_.size();i++)
        {
            distToNeighbors=0;
            currentRadius=0;
            for (int j=0;j<members_.size();j++)
            {
                float d;
                d = normScores[members_[i]->getID()][members_[j]->getID()];
                if (i!=j) {distToNeighbors += d;};
                currentRadius = (d > currentRadius) ? d : currentRadius;
            }
            if (i==0)
            {
                minDistanceSum=distToNeighbors;
                radius_=currentRadius;
            }

            if (distToNeighbors<minDistanceSum)
            {
                minDistanceSum=distToNeighbors;
                centroid_=members_[i];
                radius_=currentRadius;
            }
        }
}
