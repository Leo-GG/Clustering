/**
 * @mainpage Clustering Tools
 *
 * <b>Purpose</b>: <p>Reads a list of pairwise distances between elements
 *                 and clusters them accordingly using Hierarchical
 *                 clustering or SPICKER clustering</p>
 *
 * <b>Input</b>: <p>A list of pairwise distances</p>
 *
 * <b>Input format</b>: <p>ElementX   ElementY    distance</p>
 *
 * <b>Output</b>: <p>A list of all the clusters formed before reaching
 *                the cutoff</p>
 *
 * @author Leonardo Garma
 * @version 0.2.0 27/11/2014
 */

/**
 * @file main.cpp
 * @brief Contains the main program
 *
 * The functions for input reading, clustering and output generation
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
 * @param measureType Treat input as distances or similarities
 */
void readInput (string inpFile, int &totalNodes, vector<float> &rawScores,
                int measureType);

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

int main(int argc, char* argv[])
{

    int totalNodes;             // Number of elements to cluster
    vector<float> rawScores;    // Raw pairwise distances between elements
    int totalClusters=0;        // Counter used to assign cluster IDs

    vector< vector<float> > normScores; // Matrix of normalized distances
    vector< shared_ptr<Node> > nodeList;
    priority_queue<Link,vector<Link>,LinkComparator> linkList;
    vector<shared_ptr<Cluster> > clusterList;

    /** User input parameters **/

    string inpFile="TMData"; // Name of the file with the list of
                                    // pairwise distances
    int clusterAlg=1;    // Selected clustering algorithm
                         // Hierarchical (0) or SPICKER (1)
    int measureType=1;    // Read input as distances (0) or similarities (1)
    float cutoff=0.5908;
    int algChoice=1;

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp("-h", argv[i]))
        {
//                showHelp();
            return 0;
        }
        if (i + 1 != argc)
        { // Check that we haven't finished parsing already
            if (!strcmp("-f", argv[i]))
            {
                inpFile = argv[i + 1];
            }
            else if (!strcmp("-s", argv[i]))
            {
                clusterAlg = atoi(argv[i + 1]);
                if ( clusterAlg==0)
                {
                    algChoice=0;
                }
                else if ( clusterAlg==1)
                {
                    algChoice=1;
                }
                else
                {
                    printf ("Input Error: invalid choice of clustering algorithm\n");
                    return 1;
                }
            }
            else if (!strcmp("-m", argv[i]))
            {
                measureType = atoi(argv[i + 1]);
                if (measureType<0 || measureType>1)
                {
                    printf("Error: invalid choice of measure type\n");
                    return 1;
                }
            }
            else if (!strcmp("-d", argv[i]))
            {
                   cutoff = atof(argv[i + 1]);
            }
        }
    }

    if (measureType==1)
    {
        cutoff=(1-cutoff);
    }

    /** Clustering process **/
    readInput (inpFile, totalNodes, rawScores, measureType);

    initNodesAndClusters (totalNodes, nodeList,        // Initialize the lists
                          clusterList, totalClusters); // of Nodes and Clusters
    initScores (totalNodes,rawScores,normScores);       // Normalize the Scores

    switch (algChoice)
    {
        case 0:
            initLinks (totalNodes, normScores, // Initialize the list of Links
                       linkList, nodeList);
            doHierarchicalCutoff(linkList, clusterList, // Cluster elements
                                 totalClusters,cutoff);
            break;

        case 1:
            doSpickerCutoff(totalNodes,normScores,nodeList, clusterList,totalClusters,cutoff);
            break;
        default:
            printf ("Error: invalid choice of clustering algorithm\n");
            return 1;
    }
   // doHierarchical(linkList, clusterList, totalClusters); // Perform clustering
                                                        // using all the links
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
            if (measureType==1)
            {
                printf("radius %f ", (1-clusterList[i]->getRadius()));
            }
            else
            {
                printf("radius %f ", (clusterList[i]->getRadius()));
            }
            vector<shared_ptr<Node> > nodes = clusterList[i]->getMembers();
            printf("members %d ", nodes.size());
            if (measureType==1)
            {
                printf(", minSimilarity %f\n",(1-clusterList[i]->getMaxDistance()));
            }
            else
            {
                printf(", maxDistance %f\n",(clusterList[i]->getMaxDistance()));
            }
            printf("List of members:\n");
            for (int j = 0 ; j < nodes.size(); j++)
            {
                printf ("%d ",nodes[j]->getID());
            }
            printf("\n");
        }
    }
    printf("Total number of clusters: %d \n",activeClusters);

    return 0;
}

void readInput (string inpFile, int &totalNodes, vector<float> &rawScores, int measureType)
{

    string STRING;
	ifstream infile;
	infile.open (inpFile.c_str());
	float inputScore;
    int line=0;
        while(getline(infile,STRING))     // Read through all the lines
        {
	        vector<string> strs;          // Vector to hold words in the line
            boost::split(strs, STRING,
                         boost::is_any_of("\t+ ")); // Split the line
            if (measureType==1)
            {
                inputScore=
                ( 1-atof(strs[2].c_str()) );// Read the string with the
                                            // score as a float
            }
            else
            {
                inputScore=atof(strs[2].c_str());
            }

            rawScores.push_back(inputScore);
            line++;
        }
	infile.close();
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

void doHierarchicalCutoff(priority_queue<Link,vector<Link>,
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
                    mergeClusters(clusterList[nextLink.getNodeA()->
                                  getCluster()],
                                  clusterList[nextLink.getNodeB()->
                                  getCluster()],
                                  totalClusters++,nextLink.getDistance() );
                    clusterList.push_back(clusterC);
                }
                else
                {
                    clusterList[nextLink.getNodeA()->getCluster()]->
                    setMaxDistance(nextLink.getDistance());
                }
        }
        else
        {
            if (nextLink.getNodeA()->getCluster()==
                nextLink.getNodeB()->getCluster())
            {
                    clusterList[nextLink.getNodeA()->getCluster()]->
                    setMaxDistance(nextLink.getDistance());
            }

        }
        nextLink=linkList.top(); // Next link to check
        linkList.pop();
    }
}

void doHierarchical(priority_queue<Link,vector<Link>,
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

void doSpickerCutoff(int totalNodes, vector< vector<float> > normScores,
                     vector< shared_ptr<Node> > nodeList,
                     vector<shared_ptr<Cluster> > &clusterList,
                     int &totalClusters,float cutoff)
{
    vector< vector<float> > tempMatrix=normScores; // A copy of the distance
                                                    // matrix
    int orphans=totalNodes; // Unclustered elements
    int nextCluster=clusterList.size(); // ID for the next generated cluster

    while (orphans>0)   // Do until all elements are clustered
    {
        int maxRow=-1;  // Row on the matrix with the biggest no of neighbors
        int nbCount=0;  // Neighbor counter
        int maxNb=0;    // Maximum no of neighbors found

        vector< shared_ptr<Node> > clusterMembers;  // Elements for the next
                                                    // cluster

/* This could be improved by making a copy of the distance matrix that behaves
like a priority queue, with rows with more neighbors on the top. The second
loop, where the matrix is emptied would still be necessary. */

/** Check all rows on the matrix to find the one with more nbs*/
        for (int i =0 ; i<totalNodes;i++)
        {
            nbCount=0;
            for (int j = 0 ; j< totalNodes;j++)
            {
                if ( (tempMatrix[i][j]<cutoff) && (tempMatrix[i][j]>=0) )
                {
                    nbCount++;
                }
            }
            maxRow = nbCount >= maxNb ? i : maxRow;
            maxNb = nbCount >= maxNb ? nbCount : maxNb;
        }

/* This second loop could be improved by actually popping the matrix elements,
reducing its size and thus making the next iteration faster */

/** Push the elements of maxRow above the threshold to an array of pointers
 to Nodes and make a cluster out of them. Empty the matrix by setting the
 values to -1 */
        for (int i = 0 ; i<totalNodes ; i++)
        {
            shared_ptr<Node> node;
            if ( (tempMatrix[maxRow][i]<cutoff) && (tempMatrix[maxRow][i]>=0) )
            {
                /* Add each element below cutoff to the vector of cluster
                members and set the matrix value to -1*/
                node=nodeList[i];
                clusterMembers.push_back(node);
                clusterList[node->getCluster()]->setStatus();
                node->setCluster(nextCluster);
                tempMatrix[maxRow][i]=-1;
                orphans--;

                for (int j=0;j<totalNodes;j++)
                {
                    tempMatrix[j][i]=-1; // Reset values on all columns also!
                }
            }

        }
            /* Construct the new cluster. Then set maxDistance is set later*/
            shared_ptr<Cluster> newCluster (new Cluster(nextCluster,
                                                        clusterMembers,0));
            newCluster->calcMaxDistance(normScores);
           // newCluster->setCentroid(nodeList[maxRow]);
            clusterList.push_back(newCluster);
            nextCluster++;
    }

}

shared_ptr<Cluster> mergeClusters(shared_ptr<Cluster> A,shared_ptr<Cluster> B,
                                  int nextCluster,float maxDistance)
{
        vector<shared_ptr<Node> > clusterMembers;

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

