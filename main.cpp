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
 * <b>Usage</b>: <p>./ClustTools -f inputFile { -s algorithm | -m metric
 *                   | -d cutoff } </p>
 * @author Leonardo Garma
 * @version 0.2.0 2/12/2014
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
#include "input.h"
#include "clustering.h"

using namespace std;
using namespace std::tr1;



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
    bool hMenu=false;   // Show help?
    string inpFile="TMData"; // Name of the file with the list of
                                    // pairwise distances
    int clusterAlg=1;    // Selected clustering algorithm
                         // Hierarchical (0) or SPICKER (1)
    int measureType=1;    // Read input as distances (0) or similarities (1)
    float cutoff=0.4092;

    if (readParameters (argc, argv, hMenu, inpFile, clusterAlg,
                    measureType, cutoff) ) return 1;

    /** Clustering process **/
    readInput (inpFile, totalNodes, rawScores, measureType);

    initNodesAndClusters (totalNodes, nodeList,        // Initialize the lists
                          clusterList, totalClusters); // of Nodes and Clusters
    initScores (totalNodes,rawScores,normScores);       // Normalize the Scores

    switch (clusterAlg)
    {
        case 0:
            initLinks (totalNodes, normScores, // Initialize the list of Links
                       linkList, nodeList);
            doHierarchicalCutoff(linkList, clusterList, // Cluster elements
                                 totalClusters,cutoff);
            break;
        /*case 1:
            doHierarchical(linkList,                    // Perform clustering
                           clusterList, totalClusters); // using all the links
            break;
                                                        */
        case 1:
            doSpickerCutoff(totalNodes,normScores,nodeList, clusterList,
                            totalClusters,cutoff);
            break;
        default:
            printf ("Error: invalid choice of clustering algorithm\n");
            return 1;
    }

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
                printf(", minSimilarity %f\n",(1-clusterList[i]->
                                               getMaxDistance()));
            }
            else
            {
                printf(", maxDistance %f\n",(clusterList[i]->
                                             getMaxDistance()));
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
