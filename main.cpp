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
#include <limits>

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
    string inpFile="Identity_dist"; // Name of the file with the list of
                                    // pairwise distances
    int clusterAlg=3;    // Selected clustering algorithm
                         // Hierarchical (0) or SPICKER (1)
    int measureType=0;    // Read input as distances (0) or similarities (1)
    float cutoff=0.03;

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
        case 2:
            doKMeans(totalNodes,normScores,nodeList, clusterList,
                            totalClusters,cutoff);
            break;
        case 3:
            initLinks (totalNodes, normScores, // Initialize the list of Links
                       linkList, nodeList);
            doStrictHierarchicalCutoff(linkList, clusterList, // Cluster elements
                                 totalClusters,cutoff,normScores);
            break;
        case 4:
            initLinks (totalNodes, normScores, // Initialize the list of Links
                       linkList, nodeList);
            doUPGMA(linkList, clusterList, // Cluster elements
                                 totalClusters,cutoff,normScores);
            break;
        default:
            printf ("Error: invalid choice of clustering algorithm\n");
            return 1;


    }

    /** Output generation **/
    int activeClusters=0;
    int orphans=0;
    float maxIntra=-1;
    float minInter=std::numeric_limits<float>::max();
    float DI;
    float sumAvIntra=0;
    for (int i=0; i< clusterList.size();i++)
    {
        if (clusterList[i]->getStatus())
        {
            activeClusters++;
            clusterList[i]->calcCentroid(normScores);
            clusterList[i]->calcMean(normScores);
            clusterList[i]->calcMaxDistance(normScores);
            vector<shared_ptr<Node> > nodes = clusterList[i]->getMembers();
            maxIntra = (clusterList[i]->getMaxDistance()>maxIntra) ? clusterList[i]->getMaxDistance() : maxIntra;
            sumAvIntra+=(clusterList[i]->getAvDistance());
            if (clusterList[i]->getMembers().size() == 1) {orphans++;}

            /* Print out all cluster information */

            printf("Cluster %d : ", clusterList[i]->getID());
            printf("clustroid %d, ", clusterList[i]->getCentroid()->getID());
            printf("mean %d, ", clusterList[i]->getMean()->getID());
            printf("members %d ", nodes.size());
            if (measureType==1)
            {
                printf("radius %f ", (1-clusterList[i]->getRadius()));
                printf(", minSimilarity %f ",(1-clusterList[i]->
                                               getMaxDistance()));

                printf(", sumSimilarity %f ",(clusterList[i]->getPairs()
                                              -clusterList[i]->
                                              getDistanceSum()));


                printf(", avSimilarity %f ",((clusterList[i]->getPairs()
                                             -clusterList[i]->
                                             getDistanceSum())
                                             /clusterList[i]->getPairs()));
            }
            else
            {
                printf("radius %f ", (clusterList[i]->getRadius()));
                printf(", maxDistance %f ",(clusterList[i]->
                                             getMaxDistance()));

                printf(", sumDistance %f ",(clusterList[i]->
                                            getDistanceSum()));

                printf(", avDistance %f ",((clusterList[i]->
                                            getDistanceSum())
                                            /clusterList[i]->getPairs()));
            }
            printf(", List of members: ");
            for (int j = 0 ; j < nodes.size(); j++)
            {
                printf ("%d ",nodes[j]->getID());
            }
            printf("\n");


    /*        printf("Distance matrix: \n");
            for (int j = 0 ; j < nodes.size(); j++)
            {
                for (int k = 0 ; k < nodes.size(); k++)
                {
                    printf ("%f ",normScores[nodes[j]->getID()][nodes[k]->getID()]);
                }
                printf("\n");
            }
            printf("\n");*/
        }
    }
    printf("Total number of clusters: %d Orphans: %d ",activeClusters,orphans);
  //  printf("Total number of clusters: %d \nOrphans: %d \n",activeClusters,orphans);

    double minDist;
  //  printf ("Distance matrix for clusters: \n");
    int fakei=0;
    int fakej=0;
    float silhouetteSum=0;
    float silhouetteAv;
    float totalIntraSum=0;

    for (int i=0; i< clusterList.size()-1;i++)
    {
        if (clusterList[i]->getStatus())
        {
            fakei++;
            vector<shared_ptr<Node> > nodesi = clusterList[i]->getMembers();
            fakej=0;
      /*      for (int j=i+1; j< clusterList.size();j++)
            {
                minDist=std::numeric_limits<float>::max();
                float distSum=0;
                float avDist=0;*/

            for (int a = 0 ; a < nodesi.size(); a++)
            {
                float minAvInterDist=std::numeric_limits<float>::max();
                float distInterSum=0;
                float distIntraSum=0;
                float avInterDist;
                float avIntraDist;
                float maxIntraInter;


                if (nodesi.size()<2)
                {
                    avIntraDist=0;
                }
                else
                {
                    for (int c = 0 ; c < nodesi.size(); c++)
                    {
                        float d=normScores[nodesi[a]->getID()][nodesi[c]->getID()];
                        distIntraSum+=d;
                    }
                    if ((distIntraSum<=0) || (nodesi.size()==0))
                    {
                        avIntraDist=0;
                    }
                    else
                    {
                        avIntraDist=distIntraSum/(nodesi.size()-1);
                    }
                    totalIntraSum+=avIntraDist;
                }

             //   maxIntraSum=(distIntraSum)


                for (int j=i+1; j< clusterList.size();j++)
                {
                    if (clusterList[j]->getStatus())
                    {
                        fakej++;
                        vector<shared_ptr<Node> > nodesj = clusterList[j]->getMembers();

                        for (int b = 0 ; b < nodesj.size(); b++)
                        {
                            float d=normScores[nodesi[a]->getID()][nodesj[b]->getID()];
                         //   minDist = (d<minDist) ? d : minDist;
                            distInterSum+=d;
                        }
                        avInterDist=distInterSum/nodesj.size();
                        minAvInterDist=(avInterDist<minAvInterDist) ? avInterDist : minAvInterDist;
                    }
                //    if ( (minDist<0.9) && (minDist>0) )
                //    {
                //        printf ("%f ",minDist);
                //    }
                  //  minInter = (minDist<minInter) ? minDist : minInter;

                }


                maxIntraInter=(avIntraDist>minAvInterDist) ? avIntraDist : minAvInterDist;
                if (nodesi.size()>1)
                {
                    silhouetteSum+=(minAvInterDist-avIntraDist)/maxIntraInter;
                }

             //                           printf("%f %f %f %f\n",minAvInterDist,avIntraDist,maxIntraInter, silhouetteSum);
            }

          //  printf("\n");
        }
    }


    for (int i =0; i< nodeList.size()-1;i++)
    {
        for (int j =i+1; j<nodeList.size();j++)
        {
            int neighbors= (nodeList[i]->getCluster()==nodeList[j]->getCluster())? 1:0;
        //    printf("%d %d %d\n",i,j,neighbors);
        }


    }

    /* Print the ClusterID of every Node*/
  /*  printf("ClID=[");
    for (int i=0;i< nodeList.size()-1;i++)
    {
        printf ("%d; ",nodeList[i]->getCluster());
    }
    printf ("%d ",nodeList[nodeList.size()-1]->getCluster());
    printf("];\n");*/
    //    printf("%f %d\n",silhouetteSum,nodeList.size());

    silhouetteAv=silhouetteSum/nodeList.size();
    //DI=minInter/maxIntra;
    printf("Cutoff %f SumAvDist %f AvSil %f\n",cutoff,totalIntraSum,silhouetteAv);

    return 0;
}
