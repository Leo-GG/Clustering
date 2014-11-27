#ifndef NODE_H
#define NODE_H
/**
 * @file node.h
 * @brief Node class definition
 *
 * Defines the Node class and implements a few inline methods
 */

/**
 * @class Node
 * Represents each element in the clustering process
 * It holds a unique identifier for the node and the cluster to which it
 * belongs
 */
class Node
{
    private:

        int id_; // ID of the element
        int clusterId_; // ID of the Cluster to which it belongs

    public:

        /**
        * Constructor
        * @param id Unique identifier for the element
        * @param cluster Cluster identifier. Before the clustering process it
        *               is equal to the element identifier
        */
        Node (int id, int clusterId): id_(id),clusterId_(clusterId){};

        /**
        * Returns the cluster identifier
        * @return clusterId
        */
        int getCluster() const {return clusterId_;};

        /**
        * Returns the Node identifier
        * @return identifier
        */
        int getID() const {return id_;};

        /**
        * Assigns the element to a new cluster by changing its cluster_ value
        * @param cluster Identifier of the new cluster
        */
        void setCluster(int clusterId){clusterId_=clusterId;};

};

#endif
