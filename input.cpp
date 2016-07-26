/**
 * @file input.cpp
 * @brief Implementation of input-related functions
 *
 * Implements the functions that deal with the user input and input-file parsing:
 * readInput and readParameters
 */

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <tr1/memory>
#include <cmath>       /* sqrt */

using namespace std;
using namespace std::tr1;

void readInput (string inpFile, int &totalNodes,
                vector<float> &rawScores, int measureType)
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
                ( (1/atof(strs[2].c_str()))-1 );// Read the string with the
                ( (1-atof(strs[2].c_str())) );// Read the string with the
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

int readParameters (int argc, char* argv[], bool &hMenu,
                      string &inpFile,int &clusterAlg,int &measureType,
                      float &cutoff)
{
    /** Evaluate the parameters given by the user */
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp("-h", argv[i]))
        {
            hMenu=true;
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
            }
            else if (!strcmp("-m", argv[i]))
            {
                measureType = atoi(argv[i + 1]);
            }
            else if (!strcmp("-d", argv[i]))
            {
                   cutoff = atof(argv[i + 1]);
            }
        }
    }

    /** Evaluate options chosen */
    if (hMenu)
    {
        //showHelp();
        return 1;
    }
    if ( clusterAlg>4 && clusterAlg<0)
    {
        printf ("Input Error: invalid choice of clustering algorithm\n");
        return 1;
    }
    if (measureType<0 || measureType>1)
    {
        printf("Error: invalid choice of measure type\n");
        return 1;
    }
    if ( (measureType==1) && (clusterAlg!=2))
    {
        cutoff=1-cutoff;
        //cutoff=((1.0/cutoff)-1.0);
    }
    return 0;
}
