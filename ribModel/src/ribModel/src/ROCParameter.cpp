#include "../include/ROCParameter.h"

#include <math.h>
#include <ctime>
#include <iostream>
ROCParameter::ROCParameter()
{
    ROCParameter(1, 1, 100, 2, 1, true);
}

ROCParameter::ROCParameter(unsigned numMutationCategories, unsigned numSelectionCategories, unsigned numGenes, double sphi, unsigned _numMixtures, bool splitSer)
{
		numParam = ((splitSer) ? 41 : 40);
    Sphi = sphi;
    Sphi_proposed = sphi;
    bias_sphi = 0;
    std_sphi = 0.1;

    numAcceptForSphi = 0u;
    // proposal bias and std for phi values
    bias_phi = 0;

    // proposal bias and std for codon specific parameter
    bias_csp = 0;
    std_csp = 1;

    priorA = 0;
    priorB = 1;

    numMixtures = _numMixtures;

		initThetaK();
    currentMutationParameter.resize(numMutationCategories);
    proposedMutationParameter.resize(numMutationCategories);

    currentSelectionParameter.resize(numSelectionCategories);
    proposedSelectionParameter.resize(numSelectionCategories);

    for(unsigned i = 0; i < numMutationCategories; i++)
    {
        std::vector<double> tempMut(numParam, 0.0);
        tempMut[0] = -0.351;
        tempMut[1] = 0.378;
        tempMut[2] = 0.482;
        tempMut[3] = 0.002;
        tempMut[4] = 0.587;
        tempMut[5] = -0.420;
        tempMut[6] = 0.678;
        tempMut[7] = -0.088;
        tempMut[8] = 0.112;
        tempMut[9] = 0.308;
        tempMut[10] = 0.412;
        tempMut[11] = -0.179;
        tempMut[12] = 0.374;
        tempMut[13] = -0.594;
        tempMut[14] = 0.388;
        tempMut[15] = 0.853;
        tempMut[16] = 0.340;
        tempMut[17] = 0.463;
        tempMut[18] = 0.019;
        tempMut[19] = 0.355;
        tempMut[20] = -0.044;
        tempMut[21] = 0.292;
        tempMut[22] = 0.359;
        tempMut[23] = -0.173;
        tempMut[24] = -0.882;
        tempMut[25] = -0.699;
        tempMut[26] = -0.196;
        tempMut[27] = 0.396;
        tempMut[28] = 0.197;
        tempMut[29] = 0.023;
        tempMut[30] = 0.518;
        tempMut[31] = 0.360;
        tempMut[32] = -0.333;
        tempMut[33] = 0.286;
        tempMut[34] = 0.318;
        tempMut[35] = 0.009;
        tempMut[36] = 0.568;
        tempMut[37] = 0.197;
        tempMut[38] = 0.164;
        tempMut[39] = 0.183;


        currentMutationParameter[i] = tempMut;
        proposedMutationParameter[i] = tempMut;
    }
    for(unsigned i = 0; i < numSelectionCategories; i++)
    {
        std::vector<double> tempSel(numParam, 0.0);
        tempSel[0] = 0.683;
        tempSel[1] = 0.026;
        tempSel[2] = 0.991;
        tempSel[3] = 0.509;
        tempSel[4] = -0.216;
        tempSel[5] = -0.279;
        tempSel[6] = -0.360;
        tempSel[7] = 1.246;
        tempSel[8] = 0.544;
        tempSel[9] = 1.311;
        tempSel[10] = -0.253;
        tempSel[11] = 1.360;
        tempSel[12] = -0.146;
        tempSel[13] = 0.496;
        tempSel[14] = 0.416;
        tempSel[15] = 1.250;
        tempSel[16] = 0.697;
        tempSel[17] = 0.872;
        tempSel[18] = 0.547;
        tempSel[19] = -0.563;
        tempSel[20] = -0.442;
        tempSel[21] = 0.561;
        tempSel[22] = 0.564;
        tempSel[23] = -0.447;
        tempSel[24] = -0.009;
        tempSel[25] = 0.627;
        tempSel[26] = 2.549;
        tempSel[27] = 0.619;
        tempSel[28] = 1.829;
        tempSel[29] = 0.803;
        tempSel[30] = -0.001;
        tempSel[31] = 0.902;
        tempSel[32] = 0.747;
        tempSel[33] = -0.069;
        tempSel[34] = 0.87;
        tempSel[35] = 1.030;
        tempSel[36] = -0.051;
        tempSel[37] = 0.679;
        tempSel[38] = -0.394;
        tempSel[39] = 0.031;

        currentSelectionParameter[i] = tempSel;
        proposedSelectionParameter[i] = tempSel;
    }

    currentExpressionLevel.resize(numMixtures);
    proposedExpressionLevel.resize(numMixtures);
    for(unsigned i = 0; i < numMixtures; i++)
    {
        std::vector<double> tempExpr(numGenes, 0.0);
        currentExpressionLevel[i] = tempExpr;
        proposedExpressionLevel[i] = tempExpr;
    }

    std_phi.resize(numGenes);
    numAcceptForExpression.resize(numGenes);
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

void ROCParameter::initThetaK()
{
	for (unsigned i = 0; i < numMixtures; i++)
	{
		categories.push_back(thetaK()); //push a blank thetaK on the vector, then alter.
		categories[i].delM = 0;
		categories[i].delEta = 0;
	}

}

void ROCParameter::initAllTraces(unsigned samples, unsigned num_genes)
{
		std::cout <<"calling initExpressionTrace()\n";
    initExpressionTrace(samples, num_genes);
		std::cout <<"returned initExpressionTrace()\n";
		std::cout <<"calling initSphiTrace()\n";
    initSphiTrace(samples);
}
void ROCParameter::initExpressionTrace(unsigned samples, unsigned num_genes)
{
	/*
		std::cout <<"Starting initExpressionTrace\n";
		std::cout <<"Testing variables...\n";
		std::cout <<"numMixtures = " << numMixtures <<"\n";
		std::cout <<"samples = " << samples <<"\n";
		std::cout <<"num_genes = " << num_genes <<"\n";
		std::cout <<"expressionTrace size = " << expressionTrace.size() <<"\n";
		*/
		expressionTrace.resize(numMixtures);
    for(unsigned category = 0; category < numMixtures; category++)
    {
	//			std::cout <<"\n\n\nresizing expressionTrace\n";
        expressionTrace[category].resize(samples);
	//			std::cout <<"Going to loop through samples\n";
        for(unsigned i = 0; i < samples; i++)
        {
            expressionTrace[category][i].resize(num_genes);
            std::vector<double> tempExpr(num_genes, 0.0);
            expressionTrace[category][i] = tempExpr;
        }
    }
		/*
		std::cout <<"Final checks....\n";
		std::cout <<"expressionTrace size = " <<expressionTrace.size() <<"\n";
		std::cout <<"expressionTrace[0] size = " <<expressionTrace[0].size() <<"\n";
		std::cout <<"expressionTrace[0][0] size = " <<expressionTrace[0][0].size() <<"\n";
		std::cout <<"expressionTrace[0][0][0] = " <<expressionTrace[0][299][447] <<"\n";
		*/
}


void ROCParameter::resizeCategoryProbabilites(int samples)
{
	categoryProbabilities.resize(numMixtures);
	for (int i = 0; i < numMixtures; i++)
	{
		categoryProbabilities[i].resize(samples, 0.0);
		categoryProbabilities[i][0] = (1.0/numMixtures);
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

        unsigned* codonRange = SequenceSummary::AAindexToCodonRange(i);

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
    for(unsigned i = 0; i < genomeSize; i++)
    {
        index[i] = i;
        scuoValues[i] = calculateSCUO( genome.getGene(i) );
        expression[i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
    }
    quickSortPair(scuoValues, index, 0, genomeSize);
    quickSort(expression, 0, genomeSize);

    for(unsigned category = 0; category < numMixtures; category++)
    {
        for(unsigned j = 0; j < genomeSize; j++)
        {
            currentExpressionLevel[category][j] = expression[index[j]];
            std_phi[j] = 1;
            numAcceptForExpression[j] = 0u;
        }
    }
}
void ROCParameter::InitializeExpression(double sd_phi)
{
    int numGenes = currentExpressionLevel[1].size();
    for(unsigned category = 0; category < numMixtures; category++)
    {
        for(int i = 0; i < numGenes; i++)
        {
            currentExpressionLevel[category][i] = ROCParameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
            std_phi[i] = 1;
            numAcceptForExpression[i] = 0u;
        }
    }
}
void ROCParameter::InitializeExpression(double* expression)
{
    int numGenes = currentExpressionLevel[1].size();
    for(unsigned category = 0; category < numMixtures; category++)
    {
        for(int i = 0; i < numGenes; i++)
        {
            currentExpressionLevel[category][i] = expression[i];
            std_phi[i] = 1;
            numAcceptForExpression[i] = 0u;
        }
    }
}

void ROCParameter::getParameterForCategory(unsigned category, unsigned paramType, char aa, bool proposal, double* returnSet)
{

    std::vector<double> tempSet;
    if(paramType == ROCParameter::dM)
    {
        tempSet = proposal ? proposedMutationParameter[category] : currentMutationParameter[category];
    }
    else if(paramType == ROCParameter::dEta)
    {
        tempSet = proposal ? proposedSelectionParameter[category] : currentSelectionParameter[category];
    }
    else throw "Unkown parameter type: " + paramType;

    unsigned* aaRange = SequenceSummary::AAToCodonRange(aa, true);

    unsigned j = 0u;
    for(unsigned i = aaRange[0]; i < aaRange[1]; i++, j++)
    {
        returnSet[j] = tempSet[i];
    }
}

double ROCParameter::getExpressionPosteriorMean(unsigned samples, unsigned geneIndex, unsigned category)
{
		
    double posteriorMean = 0.0;
    unsigned traceLength = expressionTrace[0].size();
		if(samples <= traceLength)
    {
        unsigned start = traceLength - samples;
        for(unsigned i = start; i < traceLength; i++)
        {
            posteriorMean += expressionTrace[category][i][geneIndex];
        }
    }
    else{
			std::cout << std::string("ROCParameter::getExpressionPosteriorMean throws: Number of anticipated samples (") +
            std::to_string(samples) + std::string(") is greater than the length of the available trace (") + std::to_string(traceLength) + std::string(").\n");
    }
    return posteriorMean / (double)samples;
}
double ROCParameter::getSphiPosteriorMean(unsigned samples)
{
    double posteriorMean = 0.0;
    unsigned traceLength = sPhiTrace.size();
    if(samples <= traceLength)
    {
        unsigned start = traceLength - samples;

        for(unsigned i = start; i < traceLength; i++)
        {
            posteriorMean += sPhiTrace[i];
        }
    }
    else{
        throw std::string("ROCParameter::getSphiPosteriorMean throws: Number of anticipated samples (") +
            std::to_string(samples) + std::string(") is greater than the length of the available trace(") + std::to_string(traceLength) + std::string(").\n");
    }
    return posteriorMean / (double)samples;
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
    unsigned numExpressionLevels = numAcceptForExpression.size();
    for(unsigned i = 0; i < numExpressionLevels; i++)
    {
        double acceptanceLevel = (double)numAcceptForExpression[i] / (double)adaptationWidth;
        if(acceptanceLevel < 0.2)
        {
            std_phi[i] = std::max(0.01, std_phi[i] * 0.8);
        }
        if(acceptanceLevel > 0.3)
        {
            std_phi[i] = std::min(100.0, std_phi[i] * 1.2);
        }
        numAcceptForExpression[i] = 0u;
    }

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
		/*std::cout <<"Entering proposeExpressionLevels\n";
		std::cout <<"numMixtures = " <<numMixtures <<"\n";
		std::cout <<"numExpressionLevels = " <<numExpressionLevels <<"\n";
		std::cout <<"size of proposedExpressionLevel = " <<proposedExpressionLevel.size() <<"\n";
		std::cout <<"size of proposedExpressionLevel[0] = " <<proposedExpressionLevel[0].size() <<"\n";
		std::cout <<"size of currentExpressionLevel = " <<currentExpressionLevel.size() <<"\n";
		std::cout <<"size of currentExpressionLevel[0] = " <<currentExpressionLevel[0].size() <<"\n";
*/
		for(unsigned category = 0; category < numMixtures; category++)
    {
        for(unsigned i = 0; i < numExpressionLevels; i++)
        {
            // proposedExpressionLevel[i] = randLogNorm(currentExpressionLevel[i], std_phi[i]);
            // avoid adjusting probabilities for asymmetry of distribution
						proposedExpressionLevel[category][i] = std::exp( randNorm( std::log(currentExpressionLevel[category][i]) , std_phi[i]) );
        }
    }
	//	std::cout <<"Exiting proposeExpressionLevels\n";
}
void ROCParameter::proposeCodonSpecificParameter()
{
		std::cout <<"Entering proposeCodonSpecificParameter\n";
    unsigned numMutatationCategories = currentMutationParameter.size();
    for(unsigned i = 0; i < numMutatationCategories; i++)
    {
        proposedMutationParameter[i] = propose(currentMutationParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
    }
    unsigned numSelectionCategories = currentSelectionParameter.size();
    for(unsigned i = 0; i < numSelectionCategories; i++)
    {
        proposedSelectionParameter[i] = propose(currentSelectionParameter[i], ROCParameter::randNorm, bias_csp, std_csp);
    }
		std::cout <<"Exiting proposeCodonSpecificParameter\n";
}
std::vector<double> ROCParameter::propose(std::vector<double> currentParam, double (*proposal)(double a, double b), double A, double B)
{
    unsigned numParam = currentParam.size();
    std::vector<double> proposedParam(numParam, 0.0);
    for(unsigned i = 0; i < numParam; i++)
    {
        proposedParam[i] = (*proposal)(A + currentParam[i], B);
    }
    return proposedParam;
}


void ROCParameter::updateCodonSpecificParameter(unsigned category, unsigned paramType, char aa)
{
    // TODO: add acceptance counter

    unsigned* aaRange = SequenceSummary::AAToCodonRange(aa, true);

    unsigned j = 0u;
    for(unsigned i = aaRange[0]; i < aaRange[1]; i++, j++)
    {
        if(paramType == ROCParameter::dM)
        {
            currentMutationParameter[category][i] = proposedMutationParameter[category][i];
        }
        else if(paramType == ROCParameter::dEta)
        {
            currentSelectionParameter[category][i] = proposedSelectionParameter[category][i];
        }
        else throw "Unkown parameter type: " + paramType;
    }
}
/*)
void ROCParameter::updateCodonSpecificParameter(std::vector<double>& currentParamVector, std::vector<double> proposedParam,
        double (&prior)(double x, double a, double b), double logLikeImportanceRatio)
{
        // calculate log prior ratio
        double lprior = 0.0;
        for(unsigned i = 0; i < numParam; i++)
        {
            lprior += std::log( (*prior)(currentParamVector[i], priorA, priorB) ) - std::log( (*prior)(proposedParam[i], priorA, priorB) );
        }

        double logAcceptProb = logLikeImportanceRatio - lprior;

        // decide acceptance
        std::exponential_distribution<double> exponential(1);
        double propability = -exponential(generator);
        if(propability < logAcceptProb)
        {
            currentParamVector = proposedParam;
        }

*/


/*
* STATIC FUNCTIONS
*/

const unsigned ROCParameter::dM = 0;
const unsigned ROCParameter::dEta = 1;
std::default_random_engine ROCParameter::generator(time(0));

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

 unsigned ROCParameter::randMultinom(double* probabilities, unsigned groups)
 {
    // sort probabilities
    ROCParameter::quickSort(probabilities, 0, groups);
		
		std::cout << probabilities[0] <<"\n";
    // calculate cummulative sum to determine group boundaries
    double cumsum[groups];
    cumsum[0] = probabilities[0];
    for(unsigned i = 1u; i < groups; i++)
    {
        cumsum[i] = cumsum[i-1u] + probabilities[i];
    }
    // draw random number from U(0,1)
    double referenceValue = ( (double)std::rand() / (double)RAND_MAX );

    // check in which category the element falls
    unsigned returnValue = 0u;
    unsigned element = 0u;
    while(referenceValue <= cumsum[element])
    {
        returnValue = element;
        element++;
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
// qui  ck sort, sorting arrays a and b by a.
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



