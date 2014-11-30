# Hierarchical Clustering
A clustering tool based on pairwise distance measurements.
The input is a list of pairwise distances, with the follwing format:

*ElementX* *ElementY* *DistanceX-Y*

If the distances are not reciprocal ( d(X-Y) != d(Y-X) ), the program will compute
the harmonic average and use this value for the clustering. 
The program performs single-linkage clustering and it stops when it reaches a pair of
elements is larger than a given cutoff. 
The output is a list of clusters made below the cutoff that have not been merged into
a new cluster yet. For each cluster, the clustroid element, radius and maximum distance
between cluster elements are reported.


