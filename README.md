## Clustering Tools
Tools to perform clustering of elements based on pairwise distance/similarity measurements.
The input is a list of pairwise distances or similarities, with the follwing format:

*ElementX* *ElementY* *Distance/Similarity X-Y*

If the distances/similarities are not reciprocal ( d(X-Y) != d(Y-X) ), the program will compute
the harmonic average and use this value for the clustering. 
The program can perform two types of clustering:
-Single-linkage hierarchical clustering, stopping when it reaches a pair of
elements is larger than a given cutoff. 
-SPICKER clustering, where the element with the most neighbors within the cutoff is selected
iteratively as a cluster center with all its neighbors as cluster members.

The output is a list of clusters made below the cutoff that have not been merged into
a new cluster yet. For each cluster, the clustroid element, radius, maximum distance
between cluster elements and a list of members are reported.

# TO DO

*Allow for full Hierarchical clustering generation, generating a full dendogram
*Distribute the function implementations on main.cpp to separate files



