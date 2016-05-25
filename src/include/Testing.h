#ifndef Testing_H
#define Testing_H


#include "SequenceSummary.h"
#include "Gene.h"
#include "Genome.h"
#include "Utility.h"
#include "CovarianceMatrix.h"
#include "MCMCAlgorithm.h"
#include "base/Trace.h"


int testSequenceSummary();
int testGene();
int testGenome(std::string testFileDir);
int testUtility();
int testCovarianceMatrix();
int testParameter();
int testMCMCAlgorithm();
//int testTrace();

//Blank header
#endif // Testing_H
