/**
 * @file clustering.cpp
 * @brief Implementation of functions for the clustering process
 *
 * Implements all the functions needed to perform the clustering
 */

#include <tr1/memory>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <queue>
#include "node.h"
#include "cluster.h"
#include "link.h"
#include "link_comparator.h"
#include "clustering.h"
#include <iostream>
#include <fstream>


using namespace std;
using namespace std::tr1;


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
       // printf("Cluster %d has %d members\n",i,clusterList[i]->getMembers().size());
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
                if (
                    rawScores[i*totalNodes+j]==0 &&
                    rawScores[j*totalNodes+i]==0 ) score=0;
                else if (//i==j ||
                    rawScores[i*totalNodes+j]==0 ||
                    rawScores[j*totalNodes+i]==0 )
                {
                  //score=0;
                  score=(rawScores[i*totalNodes+j]+rawScores[j*totalNodes+i])/2;
                }

                else
                {
                    score = 2*
                    (rawScores[i*totalNodes+j]*rawScores[j*totalNodes+i])/
                    (rawScores[i*totalNodes+j]+rawScores[j*totalNodes+i]);
                }
                if (i==j){score=0;}
                //printf ("%d %d %f \n",i,j,score);
                row.push_back(score);
        }
       // printf ("\n");
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
                if (nextLink.getNodeA()->getCluster()!=
                    nextLink.getNodeB()->getCluster())
                    // If the linked elements are in different clusters, merge
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

void doStrictHierarchicalCutoff(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                         int totalClusters,float cutoff,
                         vector< vector<float> > normScores)
{
    Link nextLink=linkList.top(); // Next link to check
    linkList.pop();
    while (!linkList.empty())
    {
        if (nextLink.getDistance()<cutoff) // Use all the links
        {
            bool pairWise=true;
            vector<shared_ptr<Node> > nodesA =
            clusterList[nextLink.getNodeA()->getCluster()]->getMembers();
            vector<shared_ptr<Node> > nodesB =
            clusterList[nextLink.getNodeB()->getCluster()]->getMembers();

            for (int i=0; i < nodesA.size();i++)
            {
                for (int j=0; j < nodesB.size();j++)
                {
                    if (normScores[nodesA[i]->getID()][nodesB[j]->getID()]
                        >=cutoff)
                        {
                            pairWise=false;
                        }
                }
            }
            if (pairWise)
            {
                if (nextLink.getNodeA()->getCluster()!=
                    nextLink.getNodeB()->getCluster())
                    // If the linked elements are in different clusters, merge
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

void doUPGMA(priority_queue<Link,vector<Link>,
                        LinkComparator> &linkList,
                        vector<shared_ptr<Cluster> > &clusterList,
                         int totalClusters,float cutoff,
                         vector< vector<float> > normScores)
{
    Link nextLink=linkList.top(); // Next link to check
    linkList.pop();
    while (!linkList.empty())
    {
        if (nextLink.getDistance()<cutoff) // Use all the links
        {
            bool pairWise=false;
            vector<shared_ptr<Node> > nodesA =
            clusterList[nextLink.getNodeA()->getCluster()]->getMembers();
            vector<shared_ptr<Node> > nodesB =
            clusterList[nextLink.getNodeB()->getCluster()]->getMembers();
            float dSum=0;
            float avDist=0;

            for (int i=0; i < nodesA.size();i++)
            {
                for (int j=0; j < nodesB.size();j++)
                {
                    dSum+=normScores[nodesA[i]->getID()][nodesB[j]->getID()];
                }
            }
            avDist=dSum/(nodesA.size()*nodesB.size());
            if (avDist<cutoff)
            {
                pairWise=true;
            }
            if (pairWise)
            {
                if (nextLink.getNodeA()->getCluster()!=
                    nextLink.getNodeB()->getCluster())
                    // If the linked elements are in different clusters, merge
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
            /* Construct the new cluster. Then set maxDistance */
            shared_ptr<Cluster> newCluster (new Cluster(nextCluster,
                                                        clusterMembers,0));
            newCluster->calcMaxDistance(normScores);
           // newCluster->setCentroid(nodeList[maxRow]);
            clusterList.push_back(newCluster);
            nextCluster++;
    }

}


void doKMeans(int totalNodes, vector< vector<float> > normScores,
                     vector< shared_ptr<Node> > nodeList,
                     vector<shared_ptr<Cluster> > &clusterList,
                     int &totalClusters, float kMeans)
 {
    /* Initialize the k clusters */
    float dToMean;
    vector< shared_ptr<Node> > means;
    vector< vector<shared_ptr<Node> > > kMembers;
    // Make all clusters inactive
    for (int i=0; i<totalNodes ; i++)
    {
        clusterList[nodeList[i]->getCluster()]->setStatus();
    }

    // Create k new clusters at the beginning of the nodeList
    srand (time(NULL));
    for (int i =0; i<kMeans ; i++)
    {
        int newMean= rand() % totalNodes;
        bool unique=false;
        while(!unique)
        {
            unique=true;
            for (int j=0;j<i;j++)
            {
                if (newMean==means[j]->getID())
                {
                    unique=false;
                    newMean= rand() % totalNodes;
                }
            }

        }
        means.push_back(nodeList[newMean]);
        vector< shared_ptr<Node> > iMember;
        kMembers.push_back(iMember);
        iMember.push_back(means[i]);
        shared_ptr<Cluster> newCluster (new Cluster(i,iMember,0));
        clusterList.insert(clusterList.begin()+i,newCluster);
        means[i]->setCluster(i);
    }

    // Assign each element to the closest of the k new clusters
    for (int i=0; i<totalNodes ; i++)
    {
        float minDistToMean=std::numeric_limits<float>::max();
        shared_ptr<Node> closestMean=means[0];
        for (int k=0; k<kMeans; k++)
        {
            dToMean=normScores[nodeList[i]->getID()][means[k]->getID()];
            if (dToMean<minDistToMean)
            {
                minDistToMean=dToMean;
                closestMean=means[k];
              //  printf ("found a pair %d %d %f\n", i,k,dToMean);
            }
        }
        nodeList[i]->setCluster(closestMean->getCluster());
        kMembers[closestMean->getCluster()].push_back(nodeList[i]);
    }

    for (int i =0; i<kMeans ; i++)
    {
        //printf("Cluster %d mean %d\n",clusterList[i]->getID(),clusterList[i]->getMean()->getID());
        clusterList[i]->setMembers(kMembers[i]);
    }


    bool convergence=false;
    int iteration=0;
    do
    {
        convergence=true;
        /* Re-assign the mean */
        for (int i =0; i<kMeans ; i++)
        {
            int oldMean=clusterList[i]->getMean()->getID();
                    //printf("Cluster %d oldmean %d\n",clusterList[i]->getID(),oldMean);
            clusterList[i]->calcMean(normScores);
                    //printf("Cluster %d mean %d\n",clusterList[i]->getID(),clusterList[i]->getMean()->getID());

            if (oldMean!=clusterList[i]->getMean()->getID())
            {
                convergence=false;
            }
            means[i]=clusterList[i]->getMean();
            kMembers[i].clear();
        }
        //printf("iteration %d convergence %d \n", iteration, convergence);
        /* Update members based on new means */
        for (int i=0; i<totalNodes ; i++)
        {
            float minDistToMean=std::numeric_limits<float>::max();
            shared_ptr<Node> closestMean=means[0];
            for (int k=0; k<kMeans; k++)
            {
                dToMean=normScores[nodeList[i]->getID()][means[k]->getID()];
                if (dToMean<minDistToMean)
                {
                    minDistToMean=dToMean;
                    closestMean=means[k];
                  //  printf ("found a pair %d %d %f\n", i,k,dToMean);
                }
            }
            nodeList[i]->setCluster(closestMean->getCluster());
            kMembers[closestMean->getCluster()].push_back(nodeList[i]);
        }
        for (int i =0; i<kMeans ; i++)
        {
            clusterList[i]->setMembers(kMembers[i]);
        }
        iteration++;
    }while (!convergence);


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

        for (int j=0; j < B->getMembers().size();j++)
        {
            shared_ptr<Node> node=B->getMembers()[j];
            clusterMembers.push_back(node);
            node->setCluster(nextCluster);
        }

        A->setStatus();
        B->setStatus();
        shared_ptr<Cluster> C (new Cluster(nextCluster,clusterMembers,
                                           maxDistance));
        return C;
}
