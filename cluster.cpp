/**
 * @file cluster.cpp
 * @brief Implementation of methods for Cluster class
 *
 * This file contains the implementations of the calcCentroid and
 * calcMaxDistance methods.
 * calcCentroid assigns the cluster centroid and the cluster radius.
 * The centroid is defined as the element whose max distances to any
 * other member is the smallest.
 * calcMaxDistance calculates the maximum distance between any two two
 * members of the cluster.
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
                //if (i!=j) {distToNeighbors += d;};
                currentRadius = (d > currentRadius) ? d : currentRadius;
            }
            if (i==0)
            {
                minDistanceSum=currentRadius;
                radius_=currentRadius;
            }

            //if (distToNeighbors<minDistanceSum)
            if (currentRadius<minDistanceSum)
            {
                //minDistanceSum=distToNeighbors;
                minDistanceSum=currentRadius;
                centroid_=members_[i];
                radius_=currentRadius;
            }
        }
}

void Cluster::calcMaxDistance(vector< vector<float> > normScores)
{
    float maxDistance=0;
    float d;
    for (int i=0;i<members_.size();i++)
    {
        for (int j=0;j<members_.size();j++)
            {
                d=normScores[members_[i]->getID()][members_[j]->getID()];
                maxDistance = d > maxDistance ? d : maxDistance;
            }
    }
    maxDistance_=maxDistance;
}

