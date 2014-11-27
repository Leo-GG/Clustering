#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <tr1/memory>
using namespace std;
using namespace std::tr1;

/**
 * @file cluster.h
 * @brief Cluster class definition
 *
 * Defines the Cluster class and implements a few inline methods
 */

/**
 * @class Cluster
 * Represents each cluster generated in the process.
 * It holds a cluster identifier, a pointer to a vector of members
 * and the maximum distance within the cluster
 */

class Node; // Forward declaration of Node class

class Cluster
{
    private:

        int id_;                            // Unique identifier
        vector<shared_ptr<Node> > members_; /* A vector contaning all the
                                               elements in this cluster */
        float maxDistance_;                /* The maximum distance between
                                               members of this cluster */

    public:

        /**
        * Constructor
        * @param id Unique identifier
        * @param members Shared pointer to vector containing cluster members
        * @param maxDistance Maximum distance between any two members of
        *                   the cluster
        */
        Cluster(int id, vector<shared_ptr<Node> > &members,
                float maxDistance) : id_(id),members_(members),
                maxDistance_(maxDistance) {}

        /**
        * Returns the cluster identifier
        * @return id
        */
        int getID(){return id_;};

        /**
        * Returns the shared pointer to the vector of cluster members
        * @return members
        */
        vector<shared_ptr<Node> > getMembers(){return members_;};

        /**
        * Returns the maximum distance between cluster members
        * @return maxDistance
        */
        float getMaxDistance(){return maxDistance_;};

        /**
        * Updates the maximum distance within the cluster
        * @param maxDistance New maximum distance in the clsuter
        */
        void setMaxDistance(float maxDistance){maxDistance_=maxDistance;};
};

#endif
