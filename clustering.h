/**
 * @file clustering.h
 * @brief Definition of functions for the clustering process
 *
 * Defines all the functions needed to perform the clustering
 */

#ifndef CLUSTERING_H
#define CLUSTERING_H

/**
 * Generate a new Node from each element on the input file and add it to
 * the nodeList. Generate a new Cluster from each Node and assign it an ID
 * @param totalNodes Number of elements in the input/Nodes to be generated
 * @param nodeList Vector of shared pointers containing all the Nodes
 * @param clusterList Vector of shared pointers containing all the Clusters
 */
void initNodesAndClusters (int totalNodes, vector< shared_ptr<Node> > &nodeList,
                            vector<shared_ptr<Cluster> >  &clusterList,
                            int &totalClusters);

/**
 * Normalize the rawScores by computing harmonic average between reciprocal
 * scores or distances.
 * normScore=2*rawScore(A)*rawScore(B)/(rawScore(A)+rawScore(B))
 * Results in a score matrix, represented as a vector of vector of floats
 * @param totalNodes Total number of nodes, determines the size of the matrix
 * @param rawScores The distances taken from the input file
 * @param normScores Vector of vectors representing a normalized matrix of
 *                  distances between nodes
 */
void initScores(int totalNodes,vector<float> rawScores,
                vector< vector<float> > &normScores);


/**
 * Creates Links between each pair of nodes on the nodeList using the
 * normalized distance from the normScores matrix.
 * @param totalNodes Total number of nodes on the nodeList
 * @param normScores Vector of vectors representing a normalized matrix of
 *                  distances between nodes
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param nodeList Vector of shared pointers to the Nodes
 */
void initLinks (int totalNodes, vector< vector<float> > normScores,
                priority_queue<Link, vector<Link>, LinkComparator> &linkList,
                vector< shared_ptr<Node> > nodeList);


/**
 * Function for performing Hierarchical Clustering on the data set
 * Goes through all the Links in linkList and whenever two elements are not in
 * the same Cluster, merges their clusters into a new one. New clusters are
 * added to clusterList and the Cluster identifier on the Nodes are updated
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 */
void doHierarchical(priority_queue<Link,vector<Link>,LinkComparator> &linkList,
                  vector<shared_ptr<Cluster> > &clusterList, int totalClusters);

/**
 * Function for performing Hierarchical Clustering on the data set using cutoff
 * Clusters using all the Links in linkList that have a distance below a
 * given cutoff.
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 * @param cutoff A distance limit. Only links below it are considered for the
 *              clustering. The process stops when the limit is reached.
 */
void doHierarchicalCutoff(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                        int totalClusters,float cutoff);

/**
 * Function for performing Strict Hierarchical Clustering on the data set using
 * cutoff. Only merges clusters if all pairwise distances are smaller than the
 * given cutoff.
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 * @param cutoff A distance limit. Only links below it are considered for the
 *              clustering. The process stops when the limit is reached.
 */
void doStrictHierarchicalCutoff(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                        int totalClusters,float cutoff,
                        vector< vector<float> > normScores);

/**
 * Function for performing UPGMA on the data set using a given cutoff.
 * Merges clusters if the average pairwise distance is smaller than the
 * given cutoff.
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 * @param cutoff A distance limit. Only links below it are considered for the
 *              clustering. The process stops when the limit is reached.
 */
void doUPGMA(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                        int totalClusters,float cutoff,
                        vector< vector<float> > normScores);

/**
 * Joins two clusters A and B into a new Cluster C
 * that contains all the elements of A and all the elements of B
 * @param A Shared pointer to the first Cluster to join
 * @param B Shared pointer to the second Cluster to join
 * @param nextCluster Unique identifier for the Cluster to be generated
 * @param maxDistance Distance between the Node of A and the Node of B that are
 *                   connected
 * @return C Shared pointer to the newly generated Cluster
 */
shared_ptr<Cluster> mergeClusters(shared_ptr<Cluster> A,shared_ptr<Cluster> B,
                                  int nextCluster,float maxDistance);



/**
 * Function for performing clustering based on the SPICKER method
 * (Yang Z., Skolnick J., J Comput Chem. 2004 Apr 30;25(6):865-71)
 * It iteratively finds the Node with highest number of neighbors within
 * a given cutoff, forms a cluster with them and removes them from
 * the distance matrix.
 * @param totalNodes Total number of elements to cluster
 * @param normScore Vector of vector of floats represeting the distance matrix
 * @param nodeList Vector of shared pointers to the Nodes
 * @param clusterList Vector of shared pointers to the Clusters
 * @param totalClusters Total number of clusters
 * @param cutoff Distance cutoff used to perform the clustering
 */
void doSpickerCutoff(int totalNodes, vector< vector<float> > normScores,
                     vector< shared_ptr<Node> > nodeList,
                     vector<shared_ptr<Cluster> > &clusterList,
                     int &totalClusters, float cutoff);

/**
 * Function for k-means clustering. It initializes the clustering
 * with the first k elements as centroids and then repeats an assignment
 * and an update step until the assignments do not change.
 * @param totalNodes Total number of elements to cluster
 * @param normScore Vector of vector of floats represeting the distance matrix
 * @param nodeList Vector of shared pointers to the Nodes
 * @param clusterList Vector of shared pointers to the Clusters
 * @param totalClusters Total number of clusters
 * @param kMeans the k number of means used in the clustering
 */
void doKMeans(int totalNodes, vector< vector<float> > normScores,
                     vector< shared_ptr<Node> > nodeList,
                     vector<shared_ptr<Cluster> > &clusterList,
                     int &totalClusters, float kMeans);



#endif
