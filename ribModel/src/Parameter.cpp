#include "include/base/Parameter.h"
#include <sstream>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif
const std::string Parameter::allUnique = "allUnique";
const std::string Parameter::selectionShared = "selectionShared";
const std::string Parameter::mutationShared = "mutationShared";
std::default_random_engine Parameter::generator( (unsigned) std::time(NULL));

Parameter::Parameter()
{
	numParam = 0u;
	Sphi = 0.1;
	Sphi_proposed = 0.0;
	numAcceptForSphi = 0u;
	bias_sphi = 0.0;
	bias_phi = 0.0;
	numMutationCategories = 0u;
	numSelectionCategories = 0u;
	numMixtures = 0u;
	std_sphi = 0.1;
	maxGrouping = 22;
}


Parameter::Parameter(unsigned _maxGrouping)
{
	numParam = 0u;
	Sphi = 0.1;
	Sphi_proposed = 0.0;
	numAcceptForSphi = 0u;
	bias_sphi = 0.0;
	bias_phi = 0.0;
	numMutationCategories = 0u;
	numSelectionCategories = 0u;
	numMixtures = 0u;
	std_sphi = 0.1;
	maxGrouping = _maxGrouping;
}


Parameter& Parameter::operator=(const Parameter& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  numParam = rhs.numParam;


  Sphi = rhs.Sphi;
  Sphi_proposed = rhs.Sphi_proposed;
  numAcceptForSphi = rhs.numAcceptForSphi;
  categories = rhs.categories;

  // proposal bias and std for phi values
  bias_sphi = rhs.bias_sphi;
  std_sphi = rhs.std_sphi;

  // proposal bias and std for phi values
  bias_phi = rhs.bias_phi;
  std_phi = rhs.std_phi;

  currentSynthesisRateLevel = rhs.currentSynthesisRateLevel;
  proposedSynthesisRateLevel = rhs.proposedSynthesisRateLevel;
  numAcceptForSynthesisRate = rhs.numAcceptForSynthesisRate;

  numMutationCategories = rhs.numMutationCategories;
  numSelectionCategories = rhs.numSelectionCategories;


  numMixtures = rhs.numMixtures;

  mutationSelectionState = rhs.mutationSelectionState;
  selectionIsInMixture = rhs.selectionIsInMixture;
  mutationIsInMixture = rhs.mutationIsInMixture;
  maxGrouping = rhs.maxGrouping;
  groupList = rhs.groupList;
  mixtureAssignment = rhs.mixtureAssignment;
  categoryProbabilities = rhs.categoryProbabilities;
  return *this;
}

void Parameter::initParameterSet(double sphi, unsigned _numMixtures, std::vector<unsigned> geneAssignment, 
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix, bool splitSer, std::string _mutationSelectionState)
{
	// assign genes to mixture element
	unsigned numGenes = geneAssignment.size();
	mixtureAssignment.resize(numGenes, 0);
#ifndef STANDALONE
	//TODO:need to check index are correct, consecutive, and don't exceed numMixtures
	//possibly just use a set?
	for(unsigned i = 0u; i < numGenes; i++)
	{
		mixtureAssignment[i] = geneAssignment[i] - 1;
	}
#else
	for(unsigned i = 0u; i < numGenes; i++)
	{
		mixtureAssignment[i] = geneAssignment[i];
	}
#endif

	mutationSelectionState = _mutationSelectionState;
	numParam = ((splitSer) ? 40 : 41);
	numMixtures = _numMixtures;

	Sphi = sphi;
	Sphi_proposed = sphi;
	bias_sphi = 0;
	std_sphi = 0.1;

	numAcceptForSphi = 0u;

	// proposal bias and std for phi values
	bias_phi = 0;


	setNumMutationSelectionValues(_mutationSelectionState, mixtureDefinitionMatrix);
	mutationIsInMixture.resize(numMutationCategories);
	selectionIsInMixture.resize(numSelectionCategories);
	initCategoryDefinitions(_mutationSelectionState, mixtureDefinitionMatrix);

	categoryProbabilities.resize(numMixtures, 1.0/(double)numMixtures);

	currentSynthesisRateLevel.resize(numSelectionCategories);
	proposedSynthesisRateLevel.resize(numSelectionCategories);

	numAcceptForSynthesisRate.resize(numSelectionCategories);
	std_phi.resize(numSelectionCategories);

  for (unsigned i = 0u; i < numSelectionCategories; i++)
  {

    std::vector<double> tempExpr(numGenes, 0.0);
    currentSynthesisRateLevel[i] = tempExpr;
    proposedSynthesisRateLevel[i] = tempExpr;

    std::vector<unsigned> tempAccExpr(numGenes, 0u);
    numAcceptForSynthesisRate[i] = tempAccExpr;

    std::vector<double> tempStdPhi(numGenes, 1.0);
    std_phi[i] = tempStdPhi;
  }
}


// sort array interval from first (included) to last (excluded)!!
// quick sort, sorting arrays a and b by a.
// Elements in b corespond to a, a will be sorted and it will be assured that b will be sorted by a
void Parameter::quickSortPair(double a[], int b[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivotPair(a, b, first, last);
		quickSortPair(a, b, first, pivotElement);
		quickSortPair(a, b, pivotElement + 1, last);
	}
}


// sort array interval from first (included) to last (excluded)!!
void Parameter::quickSort(double a[], int first, int last)
{
	int pivotElement;

	if(first < last)
	{
		pivotElement = pivot(a, first, last);
		quickSort(a, first, pivotElement);
		quickSort(a, pivotElement + 1, last);
	}
}
int Parameter::pivot(double a[], int first, int last)
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


int Parameter::pivotPair(double a[], int b[], int first, int last)
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
void Parameter::swap(double& a, double& b)
{
	double temp = a;
	a = b;
	b = temp;
}
void Parameter::swap(int& a, int& b)
{
	int temp = a;
	a = b;
	b = temp;
}

bool Parameter::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		std::cerr <<"Error with Index\nINDEX: " << index <<"\n";
		std::cerr <<"MUST BE BETWEEN " << lowerbound << " & " << upperbound <<"\n";
	}

	return check;

}


void Parameter::writeBasicRestartFile(std::string filename)
{
	std::cout <<"Writing File\n";
	std::ofstream out;
	std::string output = "";
	std::ostringstream oss;
	unsigned i, j;
	out.open(filename.c_str());
	if (out.fail())
	{
		std::cerr <<"Could not open restart file for writing\n";
		std::exit(1);
	}


	oss << ">groupList:\n";
	for (i = 0; i < groupList.size(); i++) {
		oss << groupList[i];
		if ((i + 1) % 10 == 0) oss << "\n";
		else oss << " ";
	}
	if (i % 10 != 0) oss << "\n";
	oss <<">Sphi:\n" << Sphi <<"\n";
	oss <<">numParam:\n" << numParam <<"\n";
	oss <<">numMixtures:\n" << numMixtures <<"\n";
	oss <<">std_sphi:\n" << std_sphi <<"\n";
	//maybe clear the buffer	
	oss <<">std_phi:\n";
	for (i = 0; i < std_phi.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < std_phi[i].size(); j++)
		{
			oss << std_phi[i][j];
			if ((j + 1) % 10 == 0) oss <<"\n";
			else oss <<" ";
		}
		if (j % 10 != 0) oss <<"\n";
	}
	oss <<">categories:\n";
	for (i = 0; i < categories.size(); i++)
	{
		oss << categories[i].delM <<" " << categories[i].delEta <<"\n";	
	}
	
	oss <<">mixtureAssignment:\n";
	for (i = 0; i < mixtureAssignment.size(); i++)
	{
		oss << mixtureAssignment[i];
		if ((i + 1) % 50 == 0) oss <<"\n";
		else oss <<" ";
	}
	if (i % 50 != 0) oss <<"\n";
	oss <<">numMutationCategories:\n" << numMutationCategories <<"\n";
	oss <<">numSelectionCategories:\n" << numSelectionCategories <<"\n";

	oss <<">categoryProbabilities:\n";
	for (i = 0; i < categoryProbabilities.size(); i++)
	{
		oss << categoryProbabilities[i];
		if ((i + 1) % 10 == 0) oss <<"\n";
		else oss <<" ";
	}
	if (i % 10 != 0) oss <<"\n";


	oss <<">selectionIsInMixture:\n";
	for (i = 0; i < selectionIsInMixture.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < selectionIsInMixture[i].size(); j++)
		{
			oss << selectionIsInMixture[i][j] <<" ";
		}
		oss <<"\n";
	}

	oss <<">mutationIsInMixture:\n";
	for (i = 0; i < mutationIsInMixture.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < mutationIsInMixture[i].size(); j++)
		{
			oss << mutationIsInMixture[i][j] <<" ";
		}
		oss <<"\n";
	}

	oss <<">currentSynthesisRateLevel:\n";
	for (i = 0; i < currentSynthesisRateLevel.size(); i++)
	{
		oss <<"***\n";
		for (j = 0; j < currentSynthesisRateLevel[i].size(); j++)
		{
    	oss << currentSynthesisRateLevel[i][j];
    	if ((j + 1) % 10 == 0) oss <<"\n";
    	else oss <<" ";
		}
		if (j % 10 != 0) oss <<"\n";
	}

	std::cout <<"Done writing\n";
	output += oss.str();
	out << output;
	out.close();

}


void Parameter::initBaseValuesFromFile(std::string filename)
{
	std::ifstream input;
	input.open(filename.c_str());
	if (input.fail())
	{
		std::cerr <<"Could not open RestartFile.txt to initialzie base values\n";
		std::exit(1);
	}

	int cat = 0;
	std::vector<double> mat;
	std::string tmp, variableName;
	while (getline(input, tmp))
	{
		int flag;
		if (tmp[0] == '>') flag = 1;
		else if (input.eof()) flag = 2;
		else if (tmp[0] == '#') flag = 3;
		else flag = 4;

		if (flag == 1)
		{
			mat.clear();
			cat = 0;
			variableName = tmp.substr(1,tmp.size()-2);
		}
		else if (flag == 2)
		{
		}
		else if (flag == 3) //user comment, continue
		{
			continue;
		}
		else //store variable information
		{
			std::istringstream iss;
			if (variableName == "groupList") {
				std::string val;
				iss.str(tmp);
				while (iss >> val) {
					groupList.push_back(val);
				}
			} 
			else if (variableName == "Sphi") {iss.str(tmp); iss >> Sphi;}
			else if (variableName == "numParam") {iss.str(tmp); iss >> numParam;}	
			else if (variableName == "numMutationCategories") {iss.str(tmp); iss >> numMutationCategories;} 	
			else if (variableName == "numSelectionCategories") {iss.str(tmp); iss >> numSelectionCategories;}	
			else if (variableName == "numMixtures") {iss.str(tmp); iss >> numMixtures;}	
			else if (variableName == "mixtureAssignment") 
			{
				unsigned val;
				iss.str(tmp);
				while (iss >> val)
				{
					mixtureAssignment.push_back(val);
				}
			}	
			else if (variableName == "categories")
			{
				iss.str(tmp);
				mixtureDefinition K;
				iss >> K.delM;
				iss >> K.delEta;
				categories.push_back(K);
			} 	
			else if (variableName == "categoryProbabilities")
			{
				double val;
				iss.str(tmp);
				while (iss >> val)
				{
					categoryProbabilities.push_back(val);
				}
			} 	
			else if (variableName == "mutationIsInMixture")
			{
				if (tmp == "***")
				{
					mutationIsInMixture.resize(mutationIsInMixture.size() + 1);
					cat++;
				}
				else
				{
					unsigned val;
					iss.str(tmp);
					while (iss >> val)
					{
						mutationIsInMixture[cat - 1].push_back(val);
					}
				}
			} 	
			else if (variableName == "selectionIsInMixture")
			{
				if (tmp == "***")
				{
					selectionIsInMixture.resize(selectionIsInMixture.size() + 1);
					cat++;
				}
				else
				{
					unsigned val;
					iss.str(tmp);
					while (iss >> val)
					{
						selectionIsInMixture[cat - 1].push_back(val);
					}
				}
			} 	
			else if (variableName == "currentSynthesisRateLevel")
			{
        		if (tmp == "***")
				{
					currentSynthesisRateLevel.resize(currentSynthesisRateLevel.size() + 1);
          			cat++;
        		}
				else
				{
					double val;
          			iss.str(tmp);
					while (iss >> val)
          			{
            			currentSynthesisRateLevel[cat - 1].push_back(val);
          			}
				}
			}
			else if (variableName == "std_sphi")
			{
				iss.str(tmp);
				iss >> std_sphi;
			}
			else if (variableName == "std_phi")
			{
				if (tmp == "***")
				{
					std_phi.resize(std_phi.size() + 1);
					cat++;
				}
				iss.str(tmp);
				double val;
				while (iss >> val)
				{
					std_phi[cat - 1].push_back(val);
				}
			} 	
		}
	}

	input.close();
	
	//initialize all the default Parameter values now.
	Sphi_proposed = Sphi;
	numAcceptForSphi = 0u;
	bias_sphi = 0;
	bias_phi = 0;
	
	numAcceptForSynthesisRate.resize(numSelectionCategories);
	proposedSynthesisRateLevel.resize(numSelectionCategories);
	for (unsigned i = 0; i < numSelectionCategories; i++)
	{
		proposedSynthesisRateLevel[i] = currentSynthesisRateLevel[i];
		std::vector <double> tmp(currentSynthesisRateLevel[i].size(), 0.1);

		std::vector <unsigned> tmp2(currentSynthesisRateLevel[i].size(), 0u);
		numAcceptForSynthesisRate[i] = tmp2;
	}	
}







std::vector <double> Parameter::readPhiValues(std::string filename)
{
	std::size_t pos, pos2;
	std::ifstream currentFile;
	std::string tmpString;
	std::vector<double> RV;

	currentFile.open(filename);
	if (currentFile.fail())
	{
		std::cerr <<"Error opening file\n";
		std::exit(1);
	}

	currentFile >> tmpString; //trash the first line, no info given.


	while (currentFile >> tmpString)
	{
		pos = tmpString.find(",");
		pos2 = tmpString.find(",", pos + 1);
		if (pos != std::string::npos && pos2 != std::string::npos)
		{
			std::string val = tmpString.substr(pos + 1, pos2 - (pos + 1));
			//RV.push_back(std::stod(val));
			RV.push_back(std::atof(val.c_str()));
		}
	}

	return RV;
}


void Parameter::setNumMutationSelectionValues(std::string _mutationSelectionState, std::vector<std::vector<unsigned>> mixtureDefinitionMatrix)
{
	if (!mixtureDefinitionMatrix.empty())
	{
		//sets allow only the unique numbers to be added.
		//at the end, the size of the set is equal to the number
		//of unique categories.
		std::set<unsigned> delMCounter;
		std::set<unsigned> delEtaCounter;

		for (unsigned i = 0u; i < numMixtures; i++)
		{
			delMCounter.insert(mixtureDefinitionMatrix[i][0] - 1);
			delEtaCounter.insert(mixtureDefinitionMatrix[i][1] - 1);
		}
		numMutationCategories = delMCounter.size();
		numSelectionCategories = delEtaCounter.size();
	}
	else if (_mutationSelectionState == selectionShared)
	{
		numMutationCategories = numMixtures;
		numSelectionCategories = 1u;
	}
	else if (_mutationSelectionState == mutationShared)
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

void Parameter::initCategoryDefinitions(std::string _mutationSelectionState, std::vector<std::vector<unsigned>> mixtureDefinitionMatrix)
{
	std::set<unsigned> delMCounter;
	std::set<unsigned> delEtaCounter;

	for (unsigned i = 0; i < numMixtures; i++)
	{
		categories.push_back(mixtureDefinition()); //push a blank mixtureDefinition on the vector, then alter.
		if (!mixtureDefinitionMatrix.empty())
		{
			categories[i].delM = mixtureDefinitionMatrix[i][0] - 1;
			categories[i].delEta = mixtureDefinitionMatrix[i][1] - 1; //need check for negative and consecutive checks
			mutationIsInMixture[mixtureDefinitionMatrix[i][0] - 1].push_back(i);
			selectionIsInMixture[mixtureDefinitionMatrix[i][1] - 1].push_back(i);
		}
		else if (_mutationSelectionState == selectionShared)
		{
			categories[i].delM = i;
			categories[i].delEta = 0;
			mutationIsInMixture[i].push_back(i);
			selectionIsInMixture[0].push_back(i);
		}
		else if (_mutationSelectionState == mutationShared)
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
	//  numMutationCategories = delMCounter.size();
	//  numSelectionCategories = delEtaCounter.size();
	//sets allow only the unique numbers to be added.
	//at the end, the size of the set is equal to the number
	//of unique categories.
}



void Parameter::printMixtureDefinitionMatrix()
{
	for (unsigned i = 0u; i < numMixtures; i++)
	{
		std::cout << categories[i].delM <<"\t" << categories[i].delEta <<"\n";
	}
}

double Parameter::getSynthesisRate(unsigned geneIndex, unsigned mixtureElement, bool proposed)
{
	unsigned category = getSelectionCategory(mixtureElement);
	return (proposed ? proposedSynthesisRateLevel[category][geneIndex] : currentSynthesisRateLevel[category][geneIndex]);
}

void Parameter::setSynthesisRate(double phi, unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	currentSynthesisRateLevel[category][geneIndex] = phi;
}

double Parameter::getSynthesisRateProposalWidth(unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	return std_phi[category][geneIndex];
}

void Parameter::updateSynthesisRate(unsigned geneIndex)
{
	for(unsigned category = 0; category < numSelectionCategories; category++)
	{
		numAcceptForSynthesisRate[category][geneIndex]++;
		currentSynthesisRateLevel[category][geneIndex] = proposedSynthesisRateLevel[category][geneIndex];
	}
}

void Parameter::updateSynthesisRate(unsigned geneIndex, unsigned mixtureElement)
{
	unsigned category = getSelectionCategory(mixtureElement);
	numAcceptForSynthesisRate[category][geneIndex]++;
	currentSynthesisRateLevel[category][geneIndex] = proposedSynthesisRateLevel[category][geneIndex];
}


void Parameter::InitializeSynthesisRate(Genome& genome, double sd_phi)
{
	unsigned genomeSize = genome.getGenomeSize();
	double* scuoValues = new double[genomeSize]();
	double* expression = new double[genomeSize]();
	int* index = new int[genomeSize]();

	for(unsigned i = 0u; i < genomeSize; i++)
	{
		index[i] = i;
		scuoValues[i] = calculateSCUO( genome.getGene(i), 22 ); //This used to be maxGrouping, but RFP model will not work that way
		expression[i] = Parameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
	}
	quickSortPair(scuoValues, index, 0, genomeSize);
	quickSort(expression, 0, genomeSize);

	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned j = 0u; j < genomeSize; j++)
		{
			currentSynthesisRateLevel[category][index[j]] = expression[j];
			//std::cout << currentSynthesisRateLevel[category][j] <<"\n";
			std_phi[category][j] = 0.1;
			numAcceptForSynthesisRate[category][j] = 0u;
		}
	}

	delete [] scuoValues;
	delete [] expression;
	delete [] index;
}





void Parameter::InitializeSynthesisRate(double sd_phi)
{
	unsigned numGenes = currentSynthesisRateLevel[1].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentSynthesisRateLevel[category][i] = Parameter::randLogNorm(-(sd_phi * sd_phi) / 2, sd_phi);
			std_phi[category][i] = 0.1;
			numAcceptForSynthesisRate[category][i] = 0u;
		}
	}
}
void Parameter::InitializeSynthesisRate(std::vector<double> expression)
{
	unsigned numGenes = currentSynthesisRateLevel[0].size();
	for(unsigned category = 0u; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numGenes; i++)
		{
			currentSynthesisRateLevel[category][i] = expression[i];
			std_phi[category][i] = 0.1;
			numAcceptForSynthesisRate[category][i] = 0u;
		}
	}

}

void Parameter::proposeSphi()
{
	Sphi_proposed = std::exp(randNorm(std::log(Sphi), std_sphi));
}

void Parameter::proposeSynthesisRateLevels()
{
	unsigned numSynthesisRateLevels = currentSynthesisRateLevel[0].size();
	for(unsigned category = 0; category < numSelectionCategories; category++)
	{
		for(unsigned i = 0u; i < numSynthesisRateLevels; i++)
		{
			// avoid adjusting probabilities for asymmetry of distribution
			proposedSynthesisRateLevel[category][i] = std::exp( randNorm( std::log(currentSynthesisRateLevel[category][i]) , std_phi[category][i]) );
		}
	}
}



unsigned Parameter::getEstimatedMixtureAssignment(unsigned samples, unsigned geneIndex)
{
	unsigned rv = 0u;
	double value = -1.0;
	std::vector <double> probabilities;
	probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex);

	for (unsigned i = 0; i < probabilities.size(); i++)
	{
		if (value < probabilities[i])
		{
			value = probabilities[i];
			rv = i;
		}
	}
	return rv;
}


// calculate SCUO values according to
// Wan et al. CodonO: a new informatics method for measuring synonymous codon usage bias within and across genomes
// International Journal of General Systems, Vol. 35, No. 1, February 2006, 109–125
// http://www.tandfonline.com/doi/pdf/10.1080/03081070500502967
double Parameter::calculateSCUO(Gene& gene, unsigned maxAA)
{
	SequenceSummary seqsum = gene.getSequenceSummary();

	double totalDegenerateAACount = 0.0;
	for(unsigned i = 0; i < maxAA; i++)
	{
		std::string curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == "X" || curAA == "M" || curAA == "W") continue;
		totalDegenerateAACount += (double)seqsum.getAACountForAA(i);
	}

	double scuoValue = 0.0;
	for(unsigned i = 0; i < maxAA; i++)
	{
		std::string curAA = seqsum.AminoAcidArray[i];
		// skip amino acids with only one codon or stop codons
		if(curAA == "X" || curAA == "M" || curAA == "W") continue;
		double numDegenerateCodons = SequenceSummary::GetNumCodonsForAA(curAA);

		double aaCount = (double)seqsum.getAACountForAA(i);
		if(aaCount == 0) continue;

		std::array<unsigned, 2> codonRange = SequenceSummary::AAIndexToCodonRange(i, false);

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


void Parameter::drawIidRandomVector(unsigned draws, double mean, double sd, double (*proposal)(double a, double b), double* randomNumbers)
{
	for(unsigned i = 0u; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(mean, sd);
	}
}
void Parameter::drawIidRandomVector(unsigned draws, double r, double (*proposal)(double r), double* randomNumbers)
{
	for(unsigned i = 0u; i < draws; i++)
	{
		randomNumbers[i] = (*proposal)(r);
	}
}


double Parameter::randNorm(double mean, double sd)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rnorm(1, mean, sd);
	rv = xx[0];
#else
	std::normal_distribution<double> distribution(mean, sd);
	rv = distribution(generator);
#endif
	return rv;
}


double Parameter::randLogNorm(double m, double s)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rlnorm(1, m, s);
	rv = xx[0];
#else
	std::lognormal_distribution<double> distribution(m, s);
	rv = distribution(generator);
#endif
	return rv;
}

double Parameter::randExp(double r)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rexp(1, r);
	rv = xx[0];
#else
	std::exponential_distribution<double> distribution(r);
	rv = distribution(generator);
#endif
	return rv;
}

double Parameter::randGamma(double shape, double rate)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = rgamma(1, shape, rate = rate);
	rv = xx[0];
#else
	std::gamma_distribution<double> distribution(shape, 1 / rate);
	rv = distribution(generator);
#endif
	return rv;
}


void Parameter::randDirichlet(double* input, unsigned numElements, double* output)
{
	// draw y_i from Gamma(a_i, 1)
	// normalize y_i such that x_i = y_i / sum(y_i)

	double sumTotal = 0.0;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	for(unsigned i = 0; i < numElements; i++)
	{
		xx = rgamma(1, input[i], 1);
		output[i] = xx[0];
		sumTotal += xx[0];
	}
#else
	for(unsigned i = 0; i < numElements; i++)
	{
		std::gamma_distribution<double> distribution(input[i], 1);
		output[i] = distribution(generator);
		sumTotal += output[i];
	}
#endif
	for(unsigned i = 0; i < numElements; i++)
	{
		output[i] = output[i] / sumTotal;
	}
}


double Parameter::randUnif(double minVal, double maxVal)
{
	double rv;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, minVal, maxVal);
	rv = xx[0];
#else
	std::uniform_real_distribution<double> distribution(minVal, maxVal);
	rv = distribution(generator);
#endif
	return rv;
}

unsigned Parameter::randMultinom(double* probabilities, unsigned mixtureElements)
{
	// calculate cummulative sum to determine group boundaries
	double* cumsum = new double[mixtureElements]();
	//std::vector<double> cumsum(groups);
	cumsum[0] = probabilities[0];

	for(unsigned i = 1u; i < mixtureElements; i++)
	{
		cumsum[i] = cumsum[i-1u] + probabilities[i];
	}
	// draw random number from U(0,1)
	double referenceValue;
#ifndef STANDALONE
	RNGScope scope;
	NumericVector xx(1);
	xx = runif(1, 0, 1);
	referenceValue = xx[0];
#else
	std::uniform_real_distribution<double> distribution(0, 1);
	referenceValue = distribution(generator);
#endif
	// check in which category the element falls
	unsigned returnValue = 0u;
	for (unsigned i = 0u; i < mixtureElements; i++)
	{
		if (referenceValue <= cumsum[i])
		{
			returnValue = i;
			break;
		}
	}
	delete [] cumsum;
	return returnValue;
}

double Parameter::densityNorm(double x, double mean, double sd)
{
	const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - mean) / sd;

	return (inv_sqrt_2pi / sd) * std::exp(-0.5 * a * a);
}
double Parameter::densityLogNorm(double x, double mean, double sd)
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


//R Wrapper functions

void Parameter::setMixtureAssignmentForGene(unsigned geneIndex, unsigned value)
{
	bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
	if (check)
	{
		mixtureAssignment[geneIndex - 1] = value;
	}
}


unsigned Parameter::getMutationCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delM + 1 : 0;
}
unsigned Parameter::getSelectionCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}
unsigned Parameter::getSynthesisRateCategoryForMixture(unsigned mixtureElement)
{
	bool check = checkIndex(mixtureElement, 1, numMixtures);
	return check ? categories[mixtureElement - 1].delEta + 1 : 0;
}


unsigned Parameter::getEstimatedMixtureAssignmentForGene(unsigned samples, unsigned geneIndex)
{
	bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
	return check ? getEstimatedMixtureAssignment(samples, geneIndex - 1) + 1 : 0;
}

std::vector<double> Parameter::getEstimatedMixtureAssignmentProbabilitiesForGene(unsigned samples, unsigned geneIndex)
{
	std::vector <double> probabilities;
	bool check = checkIndex(geneIndex, 1, mixtureAssignment.size());
	if (check)
	{
		probabilities = getEstimatedMixtureAssignmentProbabilities(samples, geneIndex - 1);
	}
	return probabilities;
}



double Parameter::getSynthesisRatePosteriorMeanByMixtureElementForGene(unsigned samples, unsigned geneIndex, unsigned mixtureElement)
{
	double rv = -1.0;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignment.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures);
	if (checkGene && checkMixtureElement)
	{
		rv = getSynthesisRatePosteriorMean(samples, geneIndex - 1, mixtureElement - 1);
	}
	return rv;
}
double Parameter::getSynthesisRateVarianceByMixtureElementForGene(unsigned samples, unsigned geneIndex, unsigned mixtureElement, bool unbiased)
{
	double rv = -1.0;
	bool checkGene = checkIndex(geneIndex, 1, mixtureAssignment.size());
	bool checkMixtureElement = checkIndex(mixtureElement, 1, numMixtures);
	if (checkGene && checkMixtureElement)
	{
		rv = getSynthesisRateVariance(samples, geneIndex - 1, mixtureElement - 1, unbiased);
	}
	return rv;
}


std::vector<double> Parameter::getCurrentSynthesisRateForMixture(unsigned mixture)
{
	bool checkMixture = checkIndex(mixture, 1, numMixtures);
	unsigned exprCat = 0u;
	if(checkMixture)
	{
		exprCat = getSynthesisRateCategory(mixture - 1);
	}else{
		std::cerr << "WARNING: Mixture element " << mixture << " NOT found. Mixture element 1 is returned instead. \n";
	}
	return currentSynthesisRateLevel[exprCat];
}

void Parameter::setGroupList(std::vector <std::string> gl)
{
	groupList.clear();
	for (unsigned i = 0; i < gl.size(); i++)
	{
		if (gl[i] == "M" || gl[i] == "W" || gl[i] == "X")
		{
			std::cerr << "Warning: Amino Acid" << gl[i] << "not recognized in ROC model\n";
		}else{
			groupList.push_back(gl[i]);
		}
	}
}

std::string Parameter::getGrouping(unsigned index)
{
	return groupList[index];
}

unsigned Parameter::getGroupListSize()
{
	return groupList.size();
}
