/**
 * @mainpage Hierarchical Clustering
 *
 * <b>Purpose</b>: <p>Reads a list of pairwise distances between elements
 *                 and clusters them accordingly using Hierarchical
 *                 clustering</p>
 *
 * <b>Input</b>: <p>A list of pairwise distances</p>
 *
 * <b>Input format</b>: <p>ElementX   ElementY    distance</p>
 *
 * <b>Output</b>: <p>A list of all the clusters formed before reaching
 *                the cutoff</p>
 *
 * @author Leonardo Garma
 * @version 0.1.2 27/11/2014
 */

/**
 * @file main.cpp
 * @brief Contains the main program
 *
 * The functions for input reading reading, clustering and output generation
 * are implemented in this file.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <queue>
#include <cmath>       /* sqrt */
#include <tr1/memory>
#include "node.h"
#include "cluster.h"
#include "link.h"
#include "link_comparator.h"

using namespace std;
using namespace std::tr1;

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
 * Reads pairwise distances from a formatted file and
 * adds them to the rawScores vector. It determines how many elements are
 * taking part in the clustering.
 *
 * @param inpFile String containing the name of the input file
 * @param totalNodes Total number of nodes, determined by the lines in the
 *                  input file
 * @param rawScores Vector with the non-normalized pairwise distances
 */
void readInput (string inpFile, int &totalNodes, vector<float> &rawScores);

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
 * Goes through all the Links in linkList and whenever two elements are not in
 * the same Cluster, merges their clusters into a new one. New clusters are
 * added to clusterList and the Cluster identifier on the Nodes are updated
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 */
void doClustering(priority_queue<Link,vector<Link>,LinkComparator> &linkList,
                  vector<shared_ptr<Cluster> > &clusterList, int totalClusters);

/**
 * Clusters using all the Links in linkList that have a distance below a
 * given cutoff.
 * @param linkList Priority Queue that holds the Links in increasing order of
 *                distance
 * @param clusterList Vector of shared pointers to all the existing clusters
 * @param totalClusters Number of existing clusters
 * @param cutoff A distance limit. Only links below it are considered for the
 *              clustering. The process stops when the limit is reached.
 */
void doClusteringCutoff(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                        int totalClusters,float cutoff);


int main()
{

    int totalNodes;             // Number of elements to cluster
    vector<float> rawScores;    // Raw pairwise distances between elements
    int totalClusters=0;        // Counter used to assign cluster IDs

    vector< vector<float> > normScores; // Matrix of normalized distances
    vector< shared_ptr<Node> > nodeList;
    priority_queue<Link,vector<Link>,LinkComparator> linkList;
    vector<shared_ptr<Cluster> > clusterList;

    /** User input parameters **/
    string inpFile = "TMData";      // Name of the file with the list of pairwise
                                // distances

    float cutoff = (1-0.5908);       // Distance cutoff TM
    //float cutoff = 0.23157;          // Distance cutoff ProBis
    //float cutoff = 0.04621;       // Distance cutoff CLICK
    //float cutoff = 0.092;       // Distance cutoff Triangle
    //float cutoff = (1-0.7214);      //Dsitance cutoff SE

    /** Clustering process **/
    readInput (inpFile, totalNodes, rawScores);

    initNodesAndClusters (totalNodes, nodeList,        // Initialize the lists
                          clusterList, totalClusters); // of Nodes and Clusters

    initScores (totalNodes,rawScores,normScores);

    initLinks (totalNodes, normScores,                  // Initialize the list
               linkList, nodeList);                     // of Links

   // doClustering(linkList, clusterList, totalClusters); // Perform clustering
                                                        // using all the links

    doClusteringCutoff(linkList, clusterList,           // Perform the
                       totalClusters,cutoff);           // clustering with
                                                        // cutoff

    /** Output generation **/
    int activeClusters=0;
    for (int i=0; i< clusterList.size();i++)
    {
        if (clusterList[i]->getStatus())
        {
            activeClusters++;
            printf("Cluster %d : ", clusterList[i]->getID());
            clusterList[i]->calcCentroid(normScores);
            printf("clustroid %d, ", clusterList[i]->getCentroid()->getID());
            printf("radius %f ", (1-clusterList[i]->getRadius()));
            vector<shared_ptr<Node> > nodes = clusterList[i]->getMembers();
            printf("members %d ", nodes.size());
//            for (int j=0;j<nodes.size();j++)
//            {
//                printf("%d ", nodes[j]->getID());
//            }
            printf(", minSimilarity %f\n",(1-clusterList[i]->getMaxDistance()));
        }
    }
    printf("Total number of clusters: %d \n",activeClusters);

    return 0;
}

void readInput (string inpFile, int &totalNodes, vector<float> &rawScores)
{

    string STRING;
	ifstream infile;
	infile.open (inpFile.c_str());
	float inputScore;
    int line=0;
        while(getline(infile,STRING))     // Read through all the lines
        {
	  //      printf("%s",STRING.c_str());  // Prints our STRING
	        vector<string> strs;          // Vector to hold words in the line
            boost::split(strs, STRING,
                         boost::is_any_of("\t+ ")); // Split the line
            inputScore=( 1 -
                        atof(strs[2].c_str()) );// Read the string with the
                                                // score as a float
     //       printf (" %f\n",inputScore);
            rawScores.push_back(inputScore);
            line++;
        }
	infile.close();

	//printf("total Lines %d total Nodes %d\n",line,(int)sqrt((double)line));
    totalNodes=(int)sqrt((double)line);
}

void initNodesAndClusters (int totalNodes,
                           vector< shared_ptr<Node> > &nodeList,
                           vector<shared_ptr<Cluster> > &clusterList,
                           int &totalClusters)
{
    for (int i=0; i <totalNodes; i++)
    {
        shared_ptr<Node> node(new Node(i,i)); // Add a new element to the
                                              // list of nodes
        nodeList.push_back(node);
        vector<shared_ptr<Node> > nodes;    // Make a list of members
                                            // of a cluster
        nodes.push_back(nodeList[i]);       // Add the newly created
                                            // node the list
        shared_ptr<Cluster> cluster              // Add a new cluster
        (new Cluster (totalClusters++,nodes,0));// to the list

        clusterList.push_back(cluster);
    }
}

void initScores(int totalNodes,vector<float> rawScores,
                vector< vector<float> > &normScores)
{
    float score;
    for (int i=0; i <(totalNodes); i++)
    {
        vector <float> row;
        for (int j=0;j<totalNodes;j++)
        {
                if (i==j ||
                    rawScores[i*totalNodes+j]==0 ||
                    rawScores[j*totalNodes+i]==0 ) score=0;
                else
                {
                    score = 2*
                    (rawScores[i*totalNodes+j]*rawScores[j*totalNodes+i])/
                    (rawScores[i*totalNodes+j]+rawScores[j*totalNodes+i]);
                }
                row.push_back(score);
        }
        normScores.push_back(row);
    }
}

void initLinks (int totalNodes, vector< vector<float> > normScores,
                priority_queue<Link, vector<Link>, LinkComparator> &linkList,
                vector< shared_ptr<Node> > nodeList)
{
    for (int i=0; i <(totalNodes-1); i++)
    {
        for (int j=i+1;j<totalNodes;j++)
        {
            linkList.push(Link(nodeList[i],nodeList[j],normScores[i][j]));
        }
    }
}

void doClusteringCutoff(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                         int totalClusters,float cutoff)
{
    Link nextLink=linkList.top(); // Next link to check
    linkList.pop();
    while (!linkList.empty())
    {
        if (nextLink.getDistance()<cutoff) // Use all the links
        {
                if (nextLink.getNodeA()->getCluster()!= // If the elements of the
                    nextLink.getNodeB()->getCluster())  // link are in different
                                                        // clusters, join the
                                                        // clusters
                {
                    shared_ptr<Cluster> clusterC=
                    mergeClusters(clusterList[nextLink.getNodeA()->getCluster()],
                                  clusterList[nextLink.getNodeB()->getCluster()],
                                  totalClusters++,nextLink.getDistance() );

                    clusterList.push_back(clusterC);
                }
                else
                {
                    clusterList.at(nextLink.getNodeA()->getCluster())->setMaxDistance(nextLink.getDistance());
                }
        }
        else
        {
            if (nextLink.getNodeA()->getCluster()==nextLink.getNodeB()->getCluster())
            {
                    clusterList.at(nextLink.getNodeA()->getCluster())->setMaxDistance(nextLink.getDistance());
            }

        }
        nextLink=linkList.top(); // Next link to check
        linkList.pop();
    }
}

void doClustering(priority_queue<Link,vector<Link>,
                  LinkComparator> &linkList,
                  vector<shared_ptr<Cluster> > &clusterList,
                  int totalClusters)
{

    while (!linkList.empty()) // Use all the links
    {
            Link nextLink=linkList.top(); // Next link to check
            linkList.pop();

            if (nextLink.getNodeA()->getCluster()!= // If the elements of the
                nextLink.getNodeB()->getCluster())  // link are in different
                                                    // clusters, join the
                                                    // clusters
            {
                shared_ptr<Cluster> clusterC=
                mergeClusters(clusterList[nextLink.getNodeA()->getCluster()],
                              clusterList[nextLink.getNodeB()->getCluster()],
                              totalClusters++,nextLink.getDistance() );

                clusterList.push_back(clusterC);
            }
    }
}

shared_ptr<Cluster> mergeClusters(shared_ptr<Cluster> A,shared_ptr<Cluster> B,
                                  int nextCluster,float maxDistance)
{
        vector<shared_ptr<Node> > clusterMembers;
        vector<shared_ptr<Node> > membersA=A->getMembers();

        for (int i=0; i < A->getMembers().size();i++)
        {
            shared_ptr<Node> node=A->getMembers()[i];
            clusterMembers.push_back(node);
            node->setCluster(nextCluster);
        }

        for (int i=0; i < B->getMembers().size();i++)
        {
            shared_ptr<Node> node=B->getMembers()[i];
            clusterMembers.push_back(node);
            node->setCluster(nextCluster);
        }

        A->setStatus();
        B->setStatus();
        shared_ptr<Cluster> C (new Cluster(nextCluster,clusterMembers,
                                           maxDistance));
        return C;
}

