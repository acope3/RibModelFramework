#include "../include/ROCParameter.h"

#include <math.h>
#include <ctime>
#include <iostream>
#include <set>
#include <fstream>
ROCParameter::ROCParameter()
{
	ROCParameter(100, 2, 1, nullptr, true, "allUnique", nullptr);
}

ROCParameter::ROCParameter(unsigned numGenes, double sphi, unsigned _numMixtures, unsigned* geneAssignment, bool splitSer, std::string _mutationSelectionState,
		unsigned thetaKMatrix[][2])
{

	// assign genes to mixture element
	mixtureAssignment.resize(numGenes, 0);
	if(geneAssignment)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			mixtureAssignment[i] = geneAssignment[i];
		}
	}

    mutationSelectionState = _mutationSelectionState;
	numAcceptForMutationAndSelection.resize(22, 0u);

	numParam = ((splitSer) ? 40 : 41);
	Sphi = sphi;
	Sphi_proposed = sphi;
	bias_sphi = 0;
	std_sphi = 0.1;

	numAcceptForSphi = 0u;
	// proposal bias and std for phi values
	bias_phi = 0;

	// proposal bias and std for codon specific parameter
	bias_csp = 0;
	std_csp.resize(numParam, 0.1);

	priorA = 0;
	priorB = 1;

	numMixtures = _numMixtures;
    setNumMutationSelectionValues(_mutationSelectionState, thetaKMatrix);
    mutationIsInMixture.resize(numMutationCategories);
    selectionIsInMixture.resize(numSelectionCategories);
	initCategoryDefinitions(_mutationSelectionState, thetaKMatrix);



	currentMutationParameter.resize(numMutationCategories);
	proposedMutationParameter.resize(numMutationCategories);

	for (unsigned i = 0u; i < numMutationCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		currentMutationParameter[i] = tmp;
		proposedMutationParameter[i] = tmp;
	}


	currentSelectionParameter.resize(numSelectionCategories);
	proposedSelectionParameter.resize(numSelectionCategories);
	categoryProbabilities.resize(numMixtures, 1.0/(double)numMixtures);

	currentExpressionLevel.resize(numSelectionCategories);
	proposedExpressionLevel.resize(numSelectionCategories);

	numAcceptForExpression.resize(numSelectionCategories);
	std_phi.resize(numSelectionCategories);
	for (unsigned i = 0u; i < numSelectionCategories; i++)
	{
		std::vector<double> tmp(numParam, 0.0);
		proposedSelectionParameter[i] = tmp;
		currentSelectionParameter[i] = tmp;

		std::vector<double> tempExpr(numGenes, 0.0);
		currentExpressionLevel[i] = tempExpr;
		proposedExpressionLevel[i] = tempExpr;

		std::vector<unsigned> tempAccExpr(numGenes, 0u);
		numAcceptForExpression[i] = tempAccExpr;

		std::vector<double> tempStdPhi(numGenes, 0.1);
		std_phi[i] = tempStdPhi;
	}
}

ROCParameter::~ROCParameter()
{
	//dtor
}

ROCParameter::ROCParameter(const ROCParameter& other)
{
	numParam = other.numParam;

	Sphi = other.Sphi;
	Aphi = other.Aphi;
	Sphi_proposed = other.Sphi_proposed;
	Aphi_proposed = other.Aphi_proposed;
	sPhiTrace = other.sPhiTrace;
	aPhiTrace = other.aPhiTrace;
	numAcceptForSphi = other.numAcceptForSphi;
	categories = other.categories;

	// proposal bias and std for phi values
	bias_sphi = other.bias_sphi;
	std_sphi = other.std_sphi;

	// proposal bias and std for phi values
	bias_phi = other.bias_phi;
	std_phi = other.std_phi;

	// proposal bias and std for codon specific parameter
	bias_csp = other.bias_csp;
	std_csp = other.std_csp;

	priorA = other.priorA;
	priorB = other.priorB;

	currentExpressionLevel = other.currentExpressionLevel;
	proposedExpressionLevel = other.proposedExpressionLevel;
	expressionTrace = other.expressionTrace;
	numAcceptForExpression = other.numAcceptForExpression;

	currentMutationParameter = other.currentMutationParameter;
	proposedMutationParameter = other.proposedMutationParameter;

	currentSelectionParameter = other.currentSelectionParameter;
	proposedSelectionParameter = other.proposedSelectionParameter;
}
ROCParameter& ROCParameter::operator=(const ROCParameter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	numParam = rhs.numParam;

	Sphi = rhs.Sphi;
	Aphi = rhs.Aphi;
	Sphi_proposed = rhs.Sphi_proposed;
	Aphi_proposed = rhs.Aphi_proposed;
	sPhiTrace = rhs.sPhiTrace;
	aPhiTrace = rhs.aPhiTrace;
	numAcceptForSphi = rhs.numAcceptForSphi;
	categories = rhs.categories;

	// proposal bias and std for phi values
	bias_sphi = rhs.bias_sphi;
	std_sphi = rhs.std_sphi;

	// proposal bias and std for phi values
	bias_phi = rhs.bias_phi;
	std_phi = rhs.std_phi;

	// proposal bias and std for codon specific parameter
	bias_csp = rhs.bias_csp;
	std_csp = rhs.std_csp;

	priorA = rhs.priorA;
	priorB = rhs.priorB;

	currentExpressionLevel = rhs.currentExpressionLevel;
	proposedExpressionLevel = rhs.proposedExpressionLevel;
	expressionTrace = rhs.expressionTrace;
	numAcceptForExpression = rhs.numAcceptForExpression;

	currentMutationParameter = rhs.currentMutationParameter;
	proposedMutationParameter = rhs.proposedMutationParameter;

	currentSelectionParameter = rhs.currentSelectionParameter;
	proposedSelectionParameter = rhs.proposedSelectionParameter;
	return *this;
}


void ROCParameter::initMutationSelectionCategories(std::string files[], unsigned numCategories, unsigned paramType)
{
	unsigned i, j;
	std::size_t pos, pos2;

	std::ifstream currentFile;
	std::string tmpString;
	std::string type;

	if (paramType == ROCParameter::dM) type = "mu";
	else type = "eta";

	for (i = 0; i < numCategories; i++)
	{
		std::vector<double> temp(numParam, 0.0);

		//open the file, make sure it opens
		currentFile.open(files[i].c_str());
		if (currentFile.fail())
		{
			std::cerr <<"Error opening file " << i <<" in the file array.\n";
			std::exit(1);
		}
		currentFile >> tmpString; //trash the first line, no info given.


		j = 0;
		while (currentFile >> tmpString)
		{
			pos = tmpString.find(",");
			pos2 = tmpString.find(",", pos + 1);
			if (pos != std::string::npos && pos2 != std::string::npos)
			{
				std::string val = tmpString.substr(pos + 1, pos2 - (pos + 1));
				if (tmpString.find(type) != std::string::npos) //mu or eta was found, depending on category
				{
					temp[j] = std::stod(val);
//					temp[j] = randNorm(temp[j], 0.3);
					j++;
					if (j == numParam) break;
				}
			}
		}
		unsigned altered = 0u;
		for (j = 0; j < categories.size(); j++)
		{
			if (paramType == ROCParameter::dM && categories[j].delM == i)
			{
				currentMutationParameter[j] = temp;
				proposedMutationParameter[j] = temp;
				altered++;
			}
			else if (paramType == ROCParameter::dEta && categories[j].delEta == i)
			{
				currentSelectionParameter[j] = temp;
				proposedSelectionParameter[j] = temp;
				altered++;
			}
			if (altered == numCategories) break; //to not access indicies out of bounds.
		}
		currentFile.close();
	}
}

void ROCParameter::setNumMutationSelectionValues(std::string mutationSelectionState, unsigned thetaKMatrix[][2])
{
    if (thetaKMatrix != nullptr)
    {
        //sets allow only the unique numbers to be added.
        //at the end, the size of the set is equal to the number
        //of unique categories.
        std::set<unsigned> delMCounter;
        std::set<unsigned> delEtaCounter;

        for (unsigned i = 0u; i < numMixtures; i++)
        {
            delMCounter.insert(thetaKMatrix[i][0] - 1);
            delEtaCounter.insert(thetaKMatrix[i][1] - 1);
        }
        numMutationCategories = delMCounter.size();
        numSelectionCategories = delEtaCounter.size();
    }
    else if (mutationSelectionState == selectionShared)
    {
        numMutationCategories = numMixtures;
        numSelectionCategories = 1u;
    }
    else if (mutationSelectionState == mutationShared)
    {
        numMutationCategories = 1u;
        numSelectionCategories = numMixtures;
    }
    else //assuming the default of allUnique
    {
        numMutationCategories = numMixtures;
        numSelectionCategories = numMixtures;
    }
}

void ROCParameter::initCategoryDefinitions(std::string mutationSelectionState, unsigned thetaKMatrix[][2])
{
	std::set<unsigned> delMCounter;
	std::set<unsigned> delEtaCounter;

	for (unsigned i = 0; i < numMixtures; i++)
	{
		categories.push_back(thetaK()); //push a blank thetaK on the vector, then alter.
		if (thetaKMatrix != nullptr)
		{
			categories[i].delM = thetaKMatrix[i][0] - 1;
			categories[i].delEta = thetaKMatrix[i][1] - 1; //need check for negative and consecutive checks
			mutationIsInMixture[thetaKMatrix[i][1] - 1].push_back(i);
			selectionIsInMixture[thetaKMatrix[i][0] - 1].push_back(i);
		}
		else if (mutationSelectionState == selectionShared)
		{
			categories[i].delM = i;
			categories[i].delEta = 0;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[0].push_back(i);
		}
		else if (mutationSelectionState == mutationShared)
		{
			categories[i].delM = 0;
			categories[i].delEta = i;
			mutationIsInMixture[0].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		else //assuming the default of allUnique
		{
			categories[i].delM = i;
			categories[i].delEta = i;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[i].push_back(i);
		}
		delMCounter.insert(categories[i].delM);
		delEtaCounter.insert(categories[i].delEta);
	}
//	numMutationCategories = delMCounter.size();
//	numSelectionCategories = delEtaCounter.size();
	//sets allow only the unique numbers to be added.
	//at the end, the size of the set is equal to the number
	//of unique categories.
}

void ROCParameter::initAllTraces(unsigned samples, unsigned num_genes)
{
	initExpressionTrace(samples, num_genes);
	initMixtureAssignmentTrace(samples, num_genes);
	initSphiTrace(samples);
	initCategoryProbabilitesTrace(samples);
	initSelectionParameterTrace(samples);
	initMutationParameterTrace(samples);
}
void ROCParameter::initExpressionTrace(unsigned samples, unsigned num_genes)
{
	expressionTrace.resize(numMixtures);
	for(unsigned category = 0; category < numMixtures; category++)
	{
		expressionTrace[category].resize(samples);
		for(unsigned i = 0; i < samples; i++)
		{
			expressionTrace[category][i].resize(num_genes);
			std::vector<double> tempExpr(num_genes, 0.0);
			expressionTrace[category][i] = tempExpr;
		}
	}
}

void ROCParameter::initMutationParameterTrace(unsigned samples)
{
	mutationParameterTrace.resize(numMutationCategories);
	for (unsigned category = 0; category < numMutationCategories; category++)
	{
		mutationParameterTrace[category].resize(samples);
		for (unsigned i = 0; i < samples; i++)
		{
			mutationParameterTrace[category][i].resize(numParam);
			std::vector <double> temp(numParam, 0.0);
			mutationParameterTrace[category][i] = temp;
		}
	}
}

void ROCParameter::initSelectionParameterTrace(unsigned samples)
{
	selectionParameterTrace.resize(numSelectionCategories);
	for (unsigned category = 0; category < numSelectionCategories; category++)
	{
		selectionParameterTrace[category].resize(samples);
		for (unsigned i = 0; i < samples; i++)
		{
			selectionParameterTrace[category][i].resize(numParam);
			std::vector <double> temp(numParam, 0.0);
			selectionParameterTrace[category][i] = temp;
		}
	}
}


void ROCParameter::initMixtureAssignmentTrace(unsigned samples, unsigned num_genes)
{
	mixtureAssignmentTrace.resize(samples);
	for (unsigned i = 0u; i < samples; i ++)
	{
		mixtureAssignmentTrace[i].resize(num_genes);
	}
}


void ROCParameter::initCategoryProbabilitesTrace(int samples)
{
	categoryProbabilitiesTrace.resize(numMixtures);
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		categoryProbabilitiesTrace[i].resize(samples, 0.0);
	}

}

// calculate SCUO values according to
// Wan et al. CodonO: a new informatics method for measuring synonymous codon usage bias within and across genomes
// International Journal of General Systems, Vol. 35, No. 1, February 2006, 109–125
// http://www.tandfonline.com/doi/pdf/10.1080/03081070500502967
double ROCParameter::calculateSCUO(Gene& gene)
{
	SequenceSummary seqsum = gene.getSequenceSummary();

	double totalDegenerateAACount = 0.0;
	for(int i = 0; i < 22; i++)
	{
		char curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
		totalDegenerateAACount += (double)seqsum.getAAcountForAA(i);
	}

	double scuoValue = 0.0;
	for(int i = 0; i < 22; i++)
	{
		char curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
		double numDegenerateCodons = SequenceSummary::GetNumCodonsForAA(curAA);

		double aaCount = (double)seqsum.getAAcountForAA(i);
		if(aaCount == 0) continue;

		unsigned codonRange[2];
		SequenceSummary::AAindexToCodonRange(i, false, codonRange);

		// calculate -sum(pij log(pij))
		double aaEntropy = 0.0;
		unsigned start = codonRange[0];
		unsigned endd = codonRange[1];
		for(unsigned k = start; k < endd; k++)
		{
			int currCodonCount = seqsum.getCodonCountForCodon(k);
			if(currCodonCount == 0) continue;
			double codonProportion = (double)currCodonCount / aaCount;
			aaEntropy += codonProportion*std::log(codonProportion);
		}
		aaEntropy = -aaEntropy;
		// calculate max entropy -log(1/n_i)
		double maxEntropyForAA = -std::log(1.0 / numDegenerateCodons);
		// get normalized difference in entropy O_i
		double normalizedEntropyDiff = (maxEntropyForAA - aaEntropy) / maxEntropyForAA;

		// calculate the composition ratio F_i
		double compositionRatio = aaCount / totalDegenerateAACount;
		// SCUO is the sum(F_i * O_i) over all aa
		scuoValue += compositionRatio * normalizedEntropyDiff;
	}
	return scuoValue;
}

void ROCParameter::InitializeExpression(Genome& genome, double sd_phi)
{
	unsigned genomeSize = genome.getGenomeSize();
	double scuoValues[genomeSize];
	double expression[genomeSize];
	int index[genomeSize];
	for(unsigned i = 0u; i < genomeSize; i++)
	{
		index[i] = i;
		scuoValues[i] = calculateSCUO( genome.getGene(i) );
		expression[i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
	}
	quickSortPair(scuoValues, index, 0, genomeSize);
	quickSort(expression, 0, genomeSize);

	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned j = 0u; j < genomeSize; j++)
		{
			currentExpressionLevel[category][j] = expression[index[j]];
			std_phi[category][j] = 0.1;
			numAcceptForExpression[category][j] = 0u;
		}
	}
}
void ROCParameter::InitializeExpression(double sd_phi)
{
	unsigned numGenes = currentExpressionLevel[1].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentExpressionLevel[category][i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
			std_phi[category][i] = 0.1;
			numAcceptForExpression[category][i] = 0u;
		}
	}
}
void ROCParameter::InitializeExpression(double* expression)
{
	unsigned numGenes = currentExpressionLevel[0].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentExpressionLevel[category][i] = expression[i];
			std_phi[category][i] = 0.1;
			numAcceptForExpression[category][i] = 0u;
		}
	}

}

void ROCParameter::getParameterForCategory(unsigned category, unsigned paramType, char aa, bool proposal, double* returnSet)
{
	std::vector<double> *tempSet;
	if(paramType == ROCParameter::dM)
	{
		tempSet = (proposal ? &proposedMutationParameter[category] : &currentMutationParameter[category]);
	}
	else if(paramType == ROCParameter::dEta)
	{
		tempSet = (proposal ? &proposedSelectionParameter[category] : &currentSelectionParameter[category]);
	}
	else throw "Unkown parameter type: " + std::to_string(paramType);
	unsigned aaRange[2];
	SequenceSummary::AAToCodonRange(aa, true, aaRange);

	unsigned j = 0u;
	for(unsigned i = aaRange[0]; i < aaRange[1]; i++, j++)
	{
		if (aa =='X')
			std::cout <<"aaRange[0]: " << aaRange[0] <<" & aaRange[1]: " << aaRange[1] <<"\n";
		returnSet[j] = tempSet -> at(i);
	}
}

double ROCParameter::getMixtureAssignmentPosteriorMean(unsigned samples, unsigned geneIndex)
{

	double posteriorMean = 0.0;
	unsigned traceLength = mixtureAssignmentTrace.size();


	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getMixtureAssignmentPosteriorMean throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace (") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += mixtureAssignmentTrace[i][geneIndex];
	}

	return posteriorMean / (double)samples;
}

double ROCParameter::getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned expressionCategory)
{

	double posteriorMean = 0.0;
	unsigned traceLength = expressionTrace[0].size();


	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getExpressionPosteriorMean throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace (") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	unsigned category = 0u;
	unsigned usedSamples = 0u;
	for(unsigned i = start; i < traceLength; i++)
	{
        category = mixtureAssignmentTrace[i][geneIndex];
        category = getExpressionCategory(category);
        if(category == expressionCategory)
        {
            posteriorMean += expressionTrace[category][i][geneIndex];
            usedSamples++;
        }
	}
	// Can return NaN if gene was never in category! But that is Ok.
	return posteriorMean / (double)usedSamples;
}
double ROCParameter::getSphiPosteriorMean(unsigned samples)
{
	double posteriorMean = 0.0;
	unsigned traceLength = sPhiTrace.size();

	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getSphiPosteriorMean throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += sPhiTrace[i];
	}
	return posteriorMean / (double)samples;
}
double ROCParameter::getMutationPosteriorMean(unsigned category, unsigned samples, unsigned paramIndex)
{
	double posteriorMean = 0.0;
	unsigned traceLength = mutationParameterTrace[category].size();

	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getMutationPosteriorMean throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
    unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += mutationParameterTrace[category][i][paramIndex];
	}
	return posteriorMean / (double)samples;
}
double ROCParameter::getSelectionPosteriorMean(unsigned category, unsigned samples, unsigned paramIndex)
{
	double posteriorMean = 0.0;
	unsigned traceLength = selectionParameterTrace[category].size();

	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getSelectionPosteriorMean throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
    unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
		posteriorMean += selectionParameterTrace[category][i][paramIndex];
	}
	return posteriorMean / (double)samples;
}

double ROCParameter::getSelectionVariance(unsigned category, unsigned samples, unsigned paramIndex, bool unbiased)
{
	unsigned traceLength = selectionParameterTrace[category].size();
	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getSelectionVariance throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}

	double posteriorMean = getSelectionPosteriorMean(category, samples, paramIndex);

    double posteriorVariance = 0.0;

    unsigned start = traceLength - samples;
    double difference = 0.0;
	for(unsigned i = start; i < traceLength; i++)
	{
        difference = selectionParameterTrace[category][i][paramIndex] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}

double ROCParameter::getMutationVariance(unsigned category, unsigned samples, unsigned paramIndex, bool unbiased)
{
	unsigned traceLength = mutationParameterTrace[category].size();
	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getMutationVariance throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}

	double posteriorMean = getMutationPosteriorMean(category, samples, paramIndex);

    double posteriorVariance = 0.0;

    unsigned start = traceLength - samples;
    double difference = 0.0;
	for(unsigned i = start; i < traceLength; i++)
	{
        difference = mutationParameterTrace[category][i][paramIndex] - posteriorMean;
		posteriorVariance += difference * difference;
	}
	double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}
double ROCParameter::getExpressionVariance(unsigned samples, unsigned geneIndex, unsigned expressionCategory, bool unbiased)
{
	unsigned traceLength = expressionTrace[0].size();
	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getExpressionVariance throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace (") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}

	double posteriorMean = getExpressionPosteriorMean(samples, geneIndex, expressionCategory);

    double posteriorVariance = 0.0;
    if(!std::isnan(posteriorMean))
    {
        unsigned start = traceLength - samples;
        unsigned category = 0u;
        double difference = 0.0;
        for(unsigned i = start; i < traceLength; i++)
        {
            category = mixtureAssignmentTrace[i][geneIndex];
            category = getExpressionCategory(category);
            difference = expressionTrace[category][i][geneIndex] - posteriorMean;
            posteriorVariance += difference * difference;
        }
    }
    double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}
double ROCParameter::getSphiVariance(unsigned samples, bool unbiased)
{
	unsigned traceLength = sPhiTrace.size();
	if(samples > traceLength)
	{
		std::cout << std::string("ROCParameter::getSphiVariance throws: Number of anticipated samples (") +
			std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").") +
			std::string("Whole trace is used for posterior estimate! \n");
		samples = traceLength;
	}
	double posteriorMean = getSphiPosteriorMean(samples);

    double posteriorVariance = 0.0;

	unsigned start = traceLength - samples;
	for(unsigned i = start; i < traceLength; i++)
	{
        double difference = sPhiTrace[i] - posteriorMean;
		posteriorVariance += difference * difference;
	}
    double normalizationTerm = unbiased ? (1/((double)samples-1.0)) : (1/(double)samples);
	return normalizationTerm * posteriorVariance;
}



void ROCParameter::adaptSphiProposalWidth(unsigned adaptationWidth)
{
	double acceptanceLevel = (double)numAcceptForSphi / (double)adaptationWidth;
	if(acceptanceLevel < 0.2)
	{
		std_sphi = std::max(0.01, std_sphi * 0.8);
	}
	if(acceptanceLevel > 0.3)
	{
		std_sphi = std::min(100.0, std_sphi * 1.2);
	}
	numAcceptForSphi = 0u;
}
void ROCParameter::adaptExpressionProposalWidth(unsigned adaptationWidth)
{
	for(unsigned cat = 0u; cat < numSelectionCategories; cat ++)
	{
        unsigned numExpressionLevels = numAcceptForExpression[cat].size();
        for(unsigned i = 0; i < numExpressionLevels; i++)
        {
            double acceptanceLevel = (double)numAcceptForExpression[cat][i] / (double)adaptationWidth;
            if(acceptanceLevel < 0.2)
            {
                std_phi[cat][i] = std::max(0.01, std_phi[cat][i] * 0.8);
            }
            if(acceptanceLevel > 0.3)
            {
                std_phi[cat][i] = std::min(10.0, std_phi[cat][i] * 1.2);
            }
            numAcceptForExpression[cat][i] = 0u;
        }
	}
}
void ROCParameter::adaptCodonSpecificParameterProposalWidth(unsigned adaptationWidth)
{
    unsigned numCSPsets = numAcceptForMutationAndSelection.size();
    for(unsigned i = 0; i < numCSPsets; i++)
	{
        if(i == 21 || i == 10 || i == 18) continue;
        double acceptanceLevel = (double)numAcceptForMutationAndSelection[i] / (double)adaptationWidth;
        unsigned codonRange[2];
        SequenceSummary::AAindexToCodonRange(i, true, codonRange);
        for(unsigned k = codonRange[0]; k < codonRange[1]; k++)
        {
            if(acceptanceLevel < 0.2)
            {
                std_csp[k] = std::max(0.01, std_csp[k] * 0.8);
            }
            if(acceptanceLevel > 0.3)
            {
                std_csp[k] = std::min(10.0, std_csp[k] * 1.2);
            }
        }
		numAcceptForMutationAndSelection[i] = 0u;
	}
}

std::vector<double> ROCParameter::getExpressionTrace(unsigned geneIndex)
{
    unsigned traceLength = expressionTrace[0].size();

    std::vector<double> returnVector(traceLength, 0.0);
    for(unsigned i = 0u; i < traceLength; i++)
    {
        unsigned category = mixtureAssignmentTrace[i][geneIndex];
        returnVector[i] =  expressionTrace[category][i][geneIndex];
    }
    return returnVector;
}


// Cedric: I decided to use a normal distribution to propose Sphi and phi instead of a lognormal because:
// 1. It is a symmetric distribution and you therefore do not have to account for the unsymmetry in jump probabilities
// 2. The one log and exp operation that have to be performed per parameter are cheaper than the operations necessary to draw from a lognormal
// 3. phi has to be on a non log scale for the likelihood evaluation thus it does not help to keep phi on th elog scale all the time
// 4. the adjusment of the likelihood by the jacobian that arises from this transformation is cheap and by grouping everything in one class it takes place more or less at the same place
void ROCParameter::proposeSPhi()
{
	//Sphi_proposed = randLogNorm(Sphi, std_sphi);
	// avoid adjusting probabilities for asymmetry of distribution
	Sphi_proposed = std::exp(randNorm(std::log(Sphi), std_sphi));
}
void ROCParameter::proposeExpressionLevels()
{
	unsigned numExpressionLevels = currentExpressionLevel[0].size();
	for(unsigned category = 0; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numExpressionLevels; i++)
		{
			// avoid adjusting probabilities for asymmetry of distribution
			proposedExpressionLevel[category][i] = std::exp( randNorm( std::log(currentExpressionLevel[category][i]) , std_phi[category][i]) );
		}
	}
}
void ROCParameter::proposeCodonSpecificParameter()
{
	for(unsigned i = 0; i < numMutationCategories; i++)
	{
		proposedMutationParameter[i] = propose(currentMutationParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
	}
    for(unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedSelectionParameter[i] = propose(currentSelectionParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
    }
}
std::vector<double> ROCParameter::propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, std::vector<double> B)
{
	unsigned numParam = currentParam.size();
	std::vector<double> proposedParam(numParam, 0.0);
	for(unsigned i = 0; i < numParam; i++)
	{
		proposedParam[i] = (*proposal)(A + currentParam[i], B[i]);
	}
	return proposedParam;
}


void ROCParameter::updateCodonSpecificParameterTrace(unsigned sample, char aa)
{
    for(unsigned category = 0; category < numMutationCategories; category++)
    {
        unsigned aaRange[2];
        SequenceSummary::AAToCodonRange(aa, true, aaRange);
        for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
        {
            mutationParameterTrace[category][sample][i] = currentMutationParameter[category][i];
        }
    }
    for(unsigned category = 0; category < numSelectionCategories; category++)
    {
        unsigned aaRange[2];
        SequenceSummary::AAToCodonRange(aa, true, aaRange);
        for (unsigned i = aaRange[0]; i < aaRange[1]; i++)
        {
            selectionParameterTrace[category][sample][i] = currentSelectionParameter[category][i];
        }
    }
}

void ROCParameter::updateCodonSpecificParameter(char aa)
{
	unsigned aaRange[2];
	SequenceSummary::AAToCodonRange(aa, true, aaRange);
	unsigned aaIndex = SequenceSummary::aaToIndex.find(aa)->second;
	numAcceptForMutationAndSelection[aaIndex]++;

	for(unsigned k = 0u; k < numMutationCategories; k++)
	{
		for(unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			currentMutationParameter[k][i] = proposedMutationParameter[k][i];
		}
	}
	for(unsigned k = 0u; k < numSelectionCategories; k++)
	{
		for(unsigned i = aaRange[0]; i < aaRange[1]; i++)
		{
			currentSelectionParameter[k][i] = proposedSelectionParameter[k][i];
		}
	}
}

void ROCParameter::readPhiValues(char *filename, double temp[])
{

	int j;
	std::size_t pos, pos2;
	std::ifstream currentFile;
	std::string tmpString;
    currentFile.open(filename);
    if (currentFile.fail())
    {
      std::cerr <<"Error opening file\n";
      std::exit(1);
		}

	currentFile >> tmpString; //trash the first line, no info given.


	j = 0;
	while (currentFile >> tmpString)
	{
		pos = tmpString.find(",");
		pos2 = tmpString.find(",", pos + 1);
		if (pos != std::string::npos && pos2 != std::string::npos)
		{
			std::string val = tmpString.substr(pos + 1, pos2 - (pos + 1));
			temp[j] = std::stod(val);
			j++;
		}
	}
}

/*
 * STATIC FUNCTIONS
 */

const unsigned ROCParameter::dM = 0;
const unsigned ROCParameter::dEta = 1;
const std::string ROCParameter::allUnique = "allUnique";
const std::string ROCParameter::selectionShared = "selectionShared";
const std::string ROCParameter::mutationShared = "mutationShared";
std::default_random_engine ROCParameter::generator(time(NULL));

// TODO return array as reference and not as return
double* ROCParameter::drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b))
{
	double randomNumbers[draws];
	for(unsigned i = 0; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(mean, sd);
	}
	return randomNumbers;
}
double* ROCParameter::drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r))
{
	double randomNumbers[draws];
	for(unsigned i = 0; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(r);
	}
	return randomNumbers;
}
double ROCParameter::randNorm(double mean, double sd)
{
	std::normal_distribution<double> distribution(mean, sd);
	return distribution(generator);
}
double ROCParameter::randLogNorm(double m, double s)
{
	std::lognormal_distribution<double> distribution(m, s);
	return distribution(generator);
}
double ROCParameter::randExp(double r)
{
	std::exponential_distribution<double> distribution(r);
	return distribution(generator);
}

void ROCParameter::randDirichlet(double* input, unsigned numElements, double* output)
{
	// draw y_i from Gamma(a_i, 1)
	// normalize y_i such that x_i = y_i / sum(y_i)

	double sumTotal = 0.0;
	for(unsigned i = 0; i < numElements; i++)
	{
		std::gamma_distribution<double> distribution(input[i], 1);
		output[i] = distribution(generator);
		sumTotal += output[i];
	}
	for(unsigned i = 0; i < numElements; i++)
	{
		output[i] = output[i] / sumTotal;
	}
}

double ROCParameter::randUnif(double minVal, double maxVal)
{
	std::uniform_real_distribution<double> distribution(minVal, maxVal);
	return distribution(generator);
}

unsigned ROCParameter::randMultinom(double* probabilities, unsigned groups)
{
	// calculate cummulative sum to determine group boundaries
	double cumsum[groups];
	cumsum[0] = probabilities[0];

	for(unsigned i = 1u; i < groups; i++)
	{
		cumsum[i] = cumsum[i-1u] + probabilities[i];
	}
	// draw random number from U(0,1)
	std::uniform_real_distribution<double> distribution(0, 1);
	double referenceValue = distribution(generator);
	// check in which category the element falls
	unsigned returnValue = 0u;
	for (unsigned i = 0u; i < groups; i++)
	{
		if (referenceValue <= cumsum[i])
		{
			returnValue = i;
			break;
		}
	}
	return returnValue;
}

double ROCParameter::densityNorm(double x, double mean, double sd)
{
	const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - mean) / sd;

	return (inv_sqrt_2pi / sd) * std::exp(-0.5 * a * a);
}
double ROCParameter::densityLogNorm(double x, double mean, double sd)
{
	double returnValue = 0.0;
	// logN is only defined for x > 0 => all values less or equal to zero have probability 0
	if(x > 0.0)
	{
		const double inv_sqrt_2pi = 0.3989422804014327;
		double a = (std::log(x) - mean) / sd;
		returnValue = (inv_sqrt_2pi / (x * sd)) * std::exp(-0.5 * a * a);
	}
	return returnValue;
}

// sort array interval from first (included) to last (excluded)!!
// quick sort, sorting arrays a and b by a.
// Elements in b corespond to a, a will be sorted and it will be assured that b will be sorted by a
void ROCParameter::quickSortPair(double a[], int b[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivotPair(a, b, first, last);
		quickSortPair(a, b, first, pivotElement - 1);
		quickSortPair(a, b, pivotElement + 1, last);
	}
}
// sort array interval from first (included) to last (excluded)!!
void ROCParameter::quickSort(double a[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivot(a, first, last);
		quickSort(a, first, pivotElement - 1);
		quickSort(a, pivotElement + 1, last);
	}
}
int ROCParameter::pivot(double a[], int first, int last)
{
	int p = first;
	double pivotElement = a[first];

	for(int i = (first + 1) ; i < last ; i++)
	{
		/* If you want to sort the list in the other order, change "<=" to ">" */
		if(a[i] <= pivotElement)
		{
			p++;
			swap(a[i], a[p]);
		}
	}
	swap(a[p], a[first]);

	return p;
}
int ROCParameter::pivotPair(double a[], int b[], int first, int last)
{
	int p = first;
	double pivotElement = a[first];

	for(int i = (first + 1) ; i < last ; i++)
	{
		/* If you want to sort the list in the other order, change "<=" to ">" */
		if(a[i] <= pivotElement)
		{
			p++;
			swap(a[i], a[p]);
			swap(b[i], b[p]);
		}
	}
	swap(a[p], a[first]);
	swap(b[p], b[first]);

	return p;
}
void ROCParameter::swap(double& a, double& b)
{
	double temp = a;
	a = b;
	b = temp;
}
void ROCParameter::swap(int& a, int& b)
{
	double temp = a;
	a = b;
	b = temp;
}



