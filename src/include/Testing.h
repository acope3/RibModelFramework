#ifndef Testing_H
#define Testing_H


#include "Utility.h"
#include "SequenceSummary.h"
#include "Gene.h"
#include "Genome.h"
#include "base/Parameter.h"
#include "CovarianceMatrix.h"
#include "MCMCAlgorithm.h"
//#include "base/Trace.h"


int testUtility();
int testSequenceSummary();
int testGene();
int testGenome(std::string testFileDir);
int testParameter();
int testCovarianceMatrix();
int testRFPParameter();
int testMCMCAlgorithm();
//int testTrace();

//Blank header
#endif // Testing_H
