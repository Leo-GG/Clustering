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
 * It holds a cluster identifier, a pointer to a vector of members,
 * the maximum distance within members of the cluster, the cluster
 * centroid (member with smallest total distance to other members),
 * the cluster radius and a boolean flag that indicates whether or not
 * the cluster has been merged into a posterior cluster.
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
        shared_ptr<Node> centroid_;         // Member with smallest distance
                                            // to all other members (sum)
        float radius_;                      // largest distance from the
                                            // centroid to a member
        bool active_;

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
                maxDistance_(maxDistance)
                {
                    centroid_=members_[0];
                    radius_=0;
                    active_=true;
                }

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

        /**
        * Calculates the cluster centroid as the member whose sum
        * of distances to all other members is the smallest.
        * It assigns the value of the cluster radius (from centroid)
        * in the process
        * @param normScores A distance matrix with all the normalized
        *                    distances between members
        */
        void calcCentroid(vector< vector<float> > normScores);

        /**
        * Returns a shared pointer to the cluster centroid
        * @return centroid A shared pointer to the centroid Node
        */
        shared_ptr<Node> getCentroid(){return centroid_;};

        /**
        * Returns the cluster radius
        * @return radius The maximum distance from the centroid Node to another
        *               member of the cluster
        */
        float getRadius(){return radius_;};

        /**
        * Changes the status of the cluster, marking it as "new" or already
        * merged into a posterior one
        */
        void setStatus(){active_=!active_;};

        /**
        * Returns the active status. True if it has not yet been merged into a
        * posterior cluster, false otherwise
        * @return status Bool that indicates whether or not the cluster was
        *               formed in the "final round" of clustering
        */
        bool getStatus(){return active_;};
};

#endif
