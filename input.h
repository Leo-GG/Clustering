/**
 * @file input.h
 * @brief Definition of input-related functions
 *
 * Defines the functions that deal with the user input and input-file parsing
 */

#ifndef INPUT_H
#define INPUT_H

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
 * Reads the input parameters. And returns the corresponding choices of
 * options.
 * @param argc Number of arguments given by the user
 * @param argv Array containing the different arguments
 * @param hMenu Boolean to decide whether or not to show the help menu
 * @param inpFile String to hold the name of the input file
 * @param clusterAlg Int to hold the choice of clustering algorithm
 * @param measureType Int to hold the choice of measure
 *                   (distances/similarities)
 * @param cutoff Float to hold the value of the cutoff  for clustering
 * @return 0 if the program can continue, 1 otherwise
 */
int readParameters (int argc, char* argv[], bool &hMenu,
                      string &inpFile,int &clusterAlg,int &measureType,
                      float &cutoff);

#endif
