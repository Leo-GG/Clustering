# Hierarchical Clustering
A clustering tool based on non-reciprocal pairwise distance measurements.
The input is a list of pairwise distances, with the follwing format:

*ElementX* *ElementY* *DistanceX-Y*

The program expects to find all possible distances measured in the input (d(X-Y) and d(Y-X)).
If the distances are not reciprocal ( d(X-Y) != d(Y-X) ), the program will normalize them by
computing the harmonic average and use the average for the clustering. 

The output is not defined yet.


