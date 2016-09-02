#include "include/Genome.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif





//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//


Genome::Genome()
{
	//ctor
}


Genome& Genome::operator=(const Genome& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	genes = rhs.genes;
	simulatedGenes = rhs.simulatedGenes;
	numGenesWithPhi = rhs.numGenesWithPhi;
	//assignment operator
	return *this;
}


Genome::~Genome()
{
	//dtor
}


bool Genome::operator==(const Genome& other) const
{
	bool match = true;

 	//Do a ! operation because only the gene comparison is implemented, not the != operator.
	if (!(this->genes == other.genes)) { match = false; }
	if (!(this->simulatedGenes == other.simulatedGenes)) { match = false; }
	if (this->numGenesWithPhi != other.numGenesWithPhi) { match = false; }

	return match;
}





//----------------------------------------//
//---------- File I/O Functions ----------//
//----------------------------------------//


void Genome::readFasta(std::string filename, bool Append) // read Fasta format sequences
{
	try
	{
		if (!Append)
			clear();
		std::ifstream Fin;
		Fin.open(filename.c_str());
		if (Fin.fail())
			my_printError("ERROR: Error in Genome::readFasta: Can not open Fasta file %\n", filename);
		else
		{
			my_print("File opened\n");
			bool fastaFormat = false;
			std::string buf;
			int newLine;

			Gene tmpGene;
			std::string tempSeq = "";
			for (;;)
			{
				// read a new line in every cycle
				std::getline(Fin, buf);
				/* The new line is marked with an integer, which corresponds
					 to one of the following cases:

					 1. New sequence line, which starts with "> ".
					 2. Sequence line, which starts with any character except for ">"
					 3. End of file, detected by function eof().
				 */

				if (buf[0] == '>' ) newLine=1;
				else if (Fin.eof()) newLine=3;
				else newLine=2;

				if ( newLine == 1 )
				{ // this is a start of a new chain.
					if (!fastaFormat)
					{
						// if it is the first chain, just build a new chain object
						tmpGene.clear();
						fastaFormat = true;
					}
					else
					{
						// otherwise, need to store the old chain first, then build a new chain
						tmpGene.setSequence(tempSeq);
						addGene(tmpGene);

						tmpGene.clear();
						tempSeq = "";
					}
					tmpGene.setDescription( buf.substr(1,buf.size()-1) );
					std::size_t pos = buf.find(" ") - 1;
					tmpGene.setId( buf.substr(1,pos) );
				}

				if ( newLine == 2 && fastaFormat )
				{ // sequence line
					tempSeq.append(buf);
				}

				if ( newLine == 3 )
				{ // end of file
					my_print("EOF reached\n");
					if ( !fastaFormat )
						throw std::string("Genome::readFasta throws: ") + std::string(filename)
							  + std::string(" is not in Fasta format.");
					else
					{
						// otherwise, need to store the old chain first, then to
						// build a new chain
						tmpGene.setSequence(tempSeq);
						addGene(tmpGene);
						break;
					}
				}
			} // end while
		} // end else
		Fin.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}

void Genome::writeFasta (std::string filename, bool simulated)
{
	try {
		std::ofstream Fout;
		Fout.open(filename.c_str());
		if (Fout.fail())
			my_printError("Error in Genome::writeFasta: Can not open output Fasta file %\n", filename);
		else
		{
			unsigned sized = simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();

			for (unsigned i = 0u; i < sized; i++)
			{
				Gene *currentGene = simulated ? &simulatedGenes[i] : &genes[i];

				Fout << ">" << currentGene->getDescription() << "\n";
				for (unsigned j = 0u; j < currentGene->length(); j++)
				{
					Fout << currentGene->getNucleotideAt(j);
					if ((j + 1) % 60 == 0) Fout << std::endl;
				}
				Fout << std::endl;
			}
		} // end else
		Fout.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}

/*
 * This function reads in ONLY based on RFP_Counts > 0.
 * When this occurs, a sequence is made that is based on the number of RFP_Counts there are.
 * This is because in RFP calculations, the sequence and number of codons is technically irrelevant.
 * What is relevant is the RFP_Counts for certain codons, and an arbitrary sequence is constructed to
 * represent this.
 *
 * See also the documentation on process_sequence in SequenceSummary.cpp
 */
void Genome::readRFPFile(std::string filename)
{
	std::ifstream Fin;
	Fin.open(filename.c_str());
	if (Fin.fail())
		my_printError("Error in Genome::readRFPFile: Can not open RFP file %\n", filename);

	std::string tmp;
	//trash the first line
	if (!std::getline(Fin, tmp))
		my_printError("Error in Genome::readRFPFile: RFP file % has no header.\n", filename);
	std::string prevID = "";
	Gene tmpGene;
	bool first = true;
	std::string seq = "";

	while (std::getline(Fin, tmp))
	{
		std::size_t pos = tmp.find(",");
		std::string ID = tmp.substr(0, pos);

		if (first)
		{
			prevID = ID;
			first = false;
		}
		if (ID != prevID)
		{
			tmpGene.setId(prevID);
			tmpGene.setDescription("No description for RFP Model");
			tmpGene.setSequence(seq);
			addGene(tmpGene, false); //add to genome
			tmpGene.clear();
			seq = "";
		}
		std::size_t pos2 = tmp.find(",", pos + 1);
		std::string value = tmp.substr(pos + 1, pos2 - (pos + 1));
		unsigned counts = (unsigned)std::atoi(value.c_str());
		pos = tmp.find(",", pos2 + 1);
		std::string codon = tmp.substr(pos + 1, 3);
		for (unsigned i = 0; i < counts; i++)
			seq += codon;

		prevID = ID;
	}

	tmpGene.setId(prevID);
	tmpGene.setDescription("No description for RFP Model");
	tmpGene.setSequence(seq);

	addGene(tmpGene, false); //add to genome

	Fin.close();
}

/* Note: As the ncodons is not preserved when the RFP file is read in or RFP is otherwise processed,
 * NA is returned for the number of codons. This still preserves functionality for future
 * readRFPFile calls, as only the RFP_Counts is important.
*/
void Genome::writeRFPFile(std::string filename, bool simulated)
{
	std::ofstream Fout;
	Fout.open(filename.c_str());
	if (Fout.fail())
		my_printError("Error in Genome::writeRFPFile: Can not open output RFP file %\n", filename);
	else
	{
		Fout << "ORF,RFP_Counts,Codon_Counts,Codon\n";
		unsigned sized = simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();

		for (unsigned geneIndex = 0; geneIndex < sized; geneIndex++)
		{
			Gene *currentGene = simulated ? &simulatedGenes[geneIndex] : &genes[geneIndex];

			for (unsigned codonIndex = 0; codonIndex < 64; codonIndex++)
			{
				std::string codon = SequenceSummary::codonArray[codonIndex];

				Fout << currentGene->getId() << ",";
				Fout << currentGene->geneData.getRFPObserved(codonIndex) << ",";
				Fout << "NA," << codon << "\n";
			}
		}
	}
	Fout.close();
}


/* readPAFile (TODO: RCPP EXPOSE VIA WRAPPER)
 * Arguments: string filename, boolean to determine if we are appending to the existing genome
 * (if not set to true, will default to clearing genome data; defaults to false)
 * Read in a PA-formatted file: GeneID,Position (1-indexed),Codon,rfp_count(s) (may be multiple)
 * The positions are not necessarily in the right order.
 * There may be more than one rfp_count, and thus the header is important.
 * TODO: Wrapped by "function name" on the R-side.
 * TODO: Should we associate RFP_count with the position? Currently, it is based on order read in.
*/
void Genome::readPAFile(std::string filename, bool Append)
{
	try {
		if (!Append) clear();

		std::ifstream Fin;
		Fin.open(filename.c_str());

		if (Fin.fail())
			my_printError("Error in Genome::readPAFile: Can not open PA file %\n", filename);
		else
		{
			// Analyze the header line
			std::string tmp;
			//std::getline(Fin, tmp);

            if (!std::getline(Fin, tmp))
                my_printError("Error in Genome::readPAFile: PA file % has no header.\n", filename);

			// Ignore first 3 commas: ID, position, codon
			std::size_t pos = tmp.find(",");
			pos = tmp.find(",", pos + 1);
			pos = tmp.find(",", pos + 1);

			// While there are more commas, there are more categories of RFP counts
			unsigned numCategories = 0;
			std::size_t pos2;
			while (pos != std::string::npos)
			{
				numCategories++;
				pos2 = tmp.find(",", pos + 1);
				pos = pos2;
			}
			unsigned tableWidth = 2 + numCategories; // Table size: position + codonID + each category

			std::string prevID = "";
			Gene tmpGene;
			bool first = true;
			int position;
			//std::string seq = "";

			/*
			std::vector <unsigned> geneSizes; // Index = gene order, value = size
			unsigned geneSize = 0;

			//std::map<std::string, unsigned> genes;
			//std::map<std::string, unsigned>::iterator git;

			//First pass: For each gene name, count its size.
			while (std::getline(Fin, tmp))
			{
				pos = tmp.find(",");
				std::string ID = tmp.substr(0, pos);

				// Position integer: Ensure that if position is negative, ignore the codon.
				pos2 = tmp.find(",", pos + 1);

				// Set to 0-indexed value for convenience of vector calculation
				position = std::atoi(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str()) - 1;

				if (position > -1) // for convenience of calculation; ambiguous positions are marked by -1.
				{
					if (first)
                    {
						prevID = ID;
						first = false;
						//genes[ID] = 0;
						//geneSize = 0;
					}
					if (ID != prevID)
                    {
						geneSizes.push_back(geneSize * 3);
						geneSize = 0;

						//genes[prevID] *= 3; //multiply by three since codons
						//genes[ID] = 0; //initialize the new genes[ID]
					}
					prevID = ID;
					//genes[ID]++;
					geneSize++;
				}
			}
			geneSizes.push_back(geneSize * 3);
			//genes[prevID] *= 3; //multiply by three since codons
			Fin.close();
			 */

			// ---------------------------------------------------------------------//
			// ----- Second pass: Now for each rfp_count category, record it. ----- //
			// ----- Second pass: Treat data as a table: push back vectors of size tableWidth. ----- //
			// ---------------------------------------------------------------------//
			//Fin.open(filename.c_str());
			//std::getline(Fin, tmp); //retrash first line
			prevID = "";
			first = true;
			//unsigned geneSizesIndex = 0;
			//std::vector <std::vector <unsigned>> RFP_counts(numCategories);

			std::vector <std::vector <unsigned>> table; // Dimensions: nRows of tableWidth-sized vectors

			//Now for each line associated with a gene ID, set the string appropriately
			while (std::getline(Fin, tmp))
			{
				pos = tmp.find(",");
				std::string ID = tmp.substr(0, pos);

				// Make pos2 find next comma to find position integer
				pos2 = tmp.find(",", pos + 1);

				// Set to 0-indexed value for convenience of vector calculation
				position = std::atoi(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str()) - 1;

				// Position integer: Ensure that if position is negative, ignore the codon.
				if (position > -1) // for convenience of calculation; ambiguous positions are marked by -1.
				{
					std::vector <unsigned> tableRow(tableWidth);
					unsigned tableIndex = 0;
					tableRow[tableIndex] = (unsigned)position;

					if (first)
					{
						prevID = ID;
						first = false;
						//git = genes.find(ID);
						//seq.resize(git->second); //multiply by three since codons
					}
					if (ID != prevID)
					{
						tmpGene.setId(prevID);
						tmpGene.setDescription("No description for PANSE Model");
						//tmpGene.setSequence(seq);
						tmpGene.setPASequence(table);
						// Based on the header line, initialize the number of RFP_count categories
						// REMINDER: Modifying the sequence summary must come after setting the sequence!
						//tmpGene.initRFP_count(numCategories);

						/*
						for (unsigned i = 0; i < numCategories; i++)
							tmpGene.setRFP_count(i, RFP_counts[i]);
						*/

						addGene(tmpGene, false); //add to genome
						tmpGene.clear();
						table.clear();
						//seq = "";

						//git = genes.find(ID);
						//seq.resize(git->second); //multiply by three since codons

						/*
						for (unsigned i = 0; i < numCategories; i++)
							RFP_counts[i].clear();
						 */
					}
					// Now find codon value: Follows prior comma, guaranteed to be of size 3
					std::string codon = tmp.substr(pos2 + 1, 3);
					codon[0] = (char)std::toupper(codon[0]);
					codon[1] = (char)std::toupper(codon[1]);
					codon[2] = (char)std::toupper(codon[2]);

					tableIndex ++;
					tableRow[tableIndex] = SequenceSummary::codonToIndex(codon);
					// Note: May be an invalid Codon read in, but this is resolved when the sequence is set and processed.
					my_print("Codon is %\n", codon);

					// Skip to end rfp_count(s), if any
					pos = tmp.find(",", pos2 + 1);
					tableIndex ++;

					// While there are more commas, there are more categories of RFP counts
					while (pos != std::string::npos)
					{
						pos2 = tmp.find(",", pos + 1);
						tableRow[tableIndex] = (unsigned) std::atoi(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str());
						//unsigned RFP_count = (unsigned) std::atoi(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str());
						//RFP_counts[categoryIndex].push_back(RFP_count);
						pos = pos2;
						tableIndex++;
					}

					// Now add codon to appropriate place in sequence
					/*
					seq[position * 3] = codon[0];
					seq[position * 3 + 1] = codon[1];
					seq[position * 3 + 2] = codon[2];
					*/

					prevID = ID;
					table.push_back(tableRow);
				}
			}

            // Ensure that at least one entry was read in
            if (prevID != "")
            {
                tmpGene.setId(prevID);
                tmpGene.setDescription("No description for PANSE Model");
                //tmpGene.setSequence(seq);
                tmpGene.setPASequence(table);
                // Based on the header line, initialize the number of RFP_count categories
                // REMINDER: Modifying the sequence summary must come after setting the sequence!
                //tmpGene.initRFP_count(numCategories);

                /*
                for (unsigned i = 0; i < numCategories; i++)
                    tmpGene.setRFP_count(i, RFP_counts[i]);
                */

                addGene(tmpGene, false); //add to genome
            }
		} // end else
		Fin.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}


//TODO: Add writePANSEFile function here


/* readObservedPhiValues
 * Arguments: string filename, bool byId
 *
 * User may choose to store Phi values either by ID (if byId is true) or by index (if false).
 * In by index, each gene read in should correspond numerically to the genome.
 * Values that less than or equal to 0 or not a number will be converted to -1 and not
 * counted as a phi value.
 * NOTE: IF AN ERROR FILE IS READ IN, numGenesWithPhi IS STILL INITIALIZED WITH 0'S.
*/
void Genome::readObservedPhiValues(std::string filename, bool byId)
{
	std::ifstream input;
	std::string tmp;
	unsigned numPhi = 0;
	bool exitfunction = false;

	input.open(filename);
	if (input.fail())
		my_printError("ERROR in Genome::readObservedPhiValues: Can not open file %\n", filename);
	else
	{
		//Trash the header line
		if (!std::getline(input, tmp))
			my_printError("Error in Genome::readObservedPhiValues File: File % has no header.\n", filename);

		if (genes.size() == 0)
			my_printError("ERROR: Genome is empty, function will not execute!\n");
		else
		{
			if (byId)
			{
				//Mapping is done so the genes can be found
				std::map<std::string, Gene *> genomeMapping;
				for (unsigned i = 0; i < genes.size(); i++)
				{
					genomeMapping.insert(make_pair(genes[i].getId(), &genes[i]));
				}

                bool first = true;
                while (std::getline(input, tmp))
				{
                    std::size_t pos = tmp.find(",");
                    std::string geneID = tmp.substr(0, pos);
                    std::map<std::string, Gene *>::iterator it;
                    it = genomeMapping.find(geneID);

                    if (it == genomeMapping.end())
                        my_printError("WARNING: Gene % not found!\n", geneID);
                    else //gene is found
                    {
                        std::string val = "";
                        bool notDone = true;
						double dval;

						// Loop over each comma-separated value
                        while (notDone)
						{
							std::size_t pos2 = tmp.find(",", pos + 1);

							// End of string, end loop
							if (pos2 == std::string::npos)
							{
								val = tmp.substr(pos + 1, tmp.size() + 1);
								notDone = false;
							}
							else
							{
								val = tmp.substr(pos + 1, pos2 - (pos + 1));
								pos = pos2;
							}
							dval = std::atof(val.c_str());

							//If the value is negative, nan, or 0, we set it to -1 and throw a warning message.
							if (dval <=  0 || std::isnan(dval))
							{
								dval = -1;
								my_printError("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.\n");
							}
							it->second->observedSynthesisRateValues.push_back(dval);
						}

						// If this is the first value, initialize the size of numGenesWithPhi
						if (first)
						{
							first = false;
							numPhi = (unsigned) it->second->observedSynthesisRateValues.size();
							numGenesWithPhi.resize(numPhi, 0);
						}
						else if (it->second->observedSynthesisRateValues.size() != numPhi)
						{
                            my_printError("ERROR: Gene % has a different number of phi values given other genes: \n", geneID);
                            my_printError("Gene % has % ", geneID, it -> second -> observedSynthesisRateValues.size());
                            my_printError("while others have %. Exiting function.\n", numPhi);
							exitfunction = true;
							for (unsigned a = 0; a < getGenomeSize(); a++)
							{
								genes[a].observedSynthesisRateValues.clear();
							}
							break;
						}
                    }
                    if (exitfunction) break;
                }
				// If number of phi values match those of other given genes, execution continues normally.
				// Proceed to increment numGenesWithPhi.
                if (!exitfunction)
				{
                    for (unsigned i = 0; i < getGenomeSize(); i++)
					{
						Gene *gene = &(getGene(i));
                        if (getGene(i).observedSynthesisRateValues.size() != numPhi)
						{
							my_printError("WARNING: Gene # % (%) does not have any phi values. ", i, gene->getId());
                            my_printError("Please check your file to make sure every gene has a phi value. Filling empty genes ");
                            my_printError("with -1's for calculations.\n");
                            gene->observedSynthesisRateValues.resize(numPhi, -1);
                        }

						// Finally increment numGenesWithPhi based on stored observedSynthesisRateValues
						for (unsigned j = 0; j < numPhi; j++)
						{
							if (gene->observedSynthesisRateValues[j] != -1)
								numGenesWithPhi[j]++;
						}
                    }
                }
            }
				// By index
            else
			{
				unsigned geneIndex = 0;
				bool first = true;

				while (std::getline(input, tmp))
				{
					if (geneIndex >= genes.size())
					{
						my_printError("ERROR: GeneIndex exceeds the number of genes in the genome. Exiting function.\n");
						break;
					}

					std::size_t pos = tmp.find(",");
					std::string val = "";
					bool notDone = true;
					double dval;

					// Loop over each comma-separated value
					while (notDone)
					{
						std::size_t pos2 = tmp.find(",", pos + 1);

						// End of string, end loop
						if (pos2 == std::string::npos)
						{
							val = tmp.substr(pos + 1, tmp.size() + 1);
							notDone = false;
						}
						else
						{
							val = tmp.substr(pos + 1, pos2 - (pos + 1));
							pos = pos2;
						}
						dval = std::atof(val.c_str());

						//If the value is negative, nan, or 0, we set it to -1 and throw a warning message.
						if (dval <= 0 || std::isnan(dval))
						{
							dval = -1;
							my_printError("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.\n");
						}
						genes[geneIndex].observedSynthesisRateValues.push_back(dval);
					}

					// If this is the first value, initialize the size of numGenesWithPhi
					if (first)
					{
						first = false;
						numPhi = (unsigned) genes[geneIndex].observedSynthesisRateValues.size();
						numGenesWithPhi.resize(numPhi, 0);
					}
					else if (numPhi != genes[geneIndex].observedSynthesisRateValues.size())
					{
                        my_printError("ERROR: Gene % has a different number of phi values given other genes: \n", geneIndex);
                        my_printError("Gene % has % ", geneIndex, genes[geneIndex].observedSynthesisRateValues.size());
                        my_printError("while others have %. Exiting function.\n", numPhi);
						exitfunction = true;
						for (unsigned a = 0; a < getGenomeSize(); a++)
						{
							genes[a].observedSynthesisRateValues.clear();
						}
						break;
					}
					geneIndex++;
					if (exitfunction) break;
				}
				// If number of phi values match those of other given genes, execution continues normally.
				// Proceed to increment numGenesWithPhi.
				if (!exitfunction)
				{
					for (unsigned i = 0; i < getGenomeSize(); i++)
					{
						Gene *gene = &(getGene(i));
						if (getGene(i).observedSynthesisRateValues.size() != numPhi)
						{
                            my_printError("WARNING: Gene # % (%) does not have any phi values. ", i, gene->getId());
                            my_printError("Please check your file to make sure every gene has a phi value. Filling empty genes ");
                            my_printError("with -1's for calculations.\n");
							gene->observedSynthesisRateValues.resize(numPhi, -1);
						}

						// Finally increment numGenesWithPhi based on stored observedSynthesisRateValues
						for (unsigned j = 0; j < numPhi; j++)
						{
							if (gene->observedSynthesisRateValues[j] != -1)
								numGenesWithPhi[j]++;
						}
					}
				}
			}
		}
		input.close();
	}
}





//------------------------------------//
//---------- Gene Functions ----------//
//------------------------------------//


void Genome::addGene(const Gene& gene, bool simulated)
{
	simulated ? simulatedGenes.push_back(gene) : genes.push_back(gene);
}


std::vector <Gene> Genome::getGenes(bool simulated)
{
	return !simulated ? genes : simulatedGenes;
}


unsigned Genome::getNumGenesWithPhiForIndex(unsigned index)
{
	return numGenesWithPhi[index];
}


Gene& Genome::getGene(unsigned index, bool simulated)
{
	return simulated ? simulatedGenes[index] : genes[index];
}


Gene& Genome::getGene(std::string id, bool simulated)
{
	Gene tempGene;
	unsigned geneIndex;

	for (geneIndex = 0; geneIndex < getGenomeSize(); geneIndex++)
	{
		tempGene = simulated ? simulatedGenes[geneIndex] : genes[geneIndex];

		if (tempGene.getId().compare(id) == 0) break;
	}
	return simulated ? simulatedGenes[geneIndex] : genes[geneIndex];
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


unsigned Genome::getGenomeSize(bool simulated)
{
	return simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();
}


void Genome::clear()
{
	genes.clear();
	simulatedGenes.clear();
	numGenesWithPhi.clear();
}


Genome Genome::getGenomeForGeneIndices(std::vector <unsigned> indices, bool simulated)
{
	Genome genome;

	for (unsigned i = 0; i < indices.size(); i++)
	{
		if (indices[i] > getGenomeSize(simulated))
		{
			my_printError("Error in Genome::getGenomeForGeneIndices. An index specified is out of bounds for the genome!\n");
			my_printError("The index % is greater than the size of the genome (%).\n", indices[i], getGenomeSize());
			my_printError("Returning empty Genome.\n");
			genome.clear();
			return genome;
		}
		else
		{
			simulated ? genome.addGene(simulatedGenes[indices[i]], true) : genome.addGene(genes[indices[i]], false);
		}
	}

	return genome;
}


std::vector<unsigned> Genome::getCodonCountsPerGene(std::string codon)
{
	std::vector<unsigned> codonCounts(genes.size());
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	for (unsigned i = 0u; i < genes.size(); i++)
	{
		Gene gene = genes[i];
		SequenceSummary *seqsum = gene.getSequenceSummary();
		codonCounts[i] = seqsum -> getCodonCountForCodon(codonIndex);
	}
	return codonCounts;
}





//---------------------------------------//
//---------- Testing Functions ----------//
//---------------------------------------//


std::vector <unsigned> Genome::getNumGenesWithPhi()
{
	return numGenesWithPhi;
}


void Genome::setNumGenesWithPhi(std::vector <unsigned> newVector)
{
	this->numGenesWithPhi = newVector;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE


bool Genome::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		my_printError("ERROR: Index % is out of bounds. Index must be between % & %\n", index, lowerbound, upperbound);
	}
	return check;
}

//NOTE: This function does the check and performs the function itself because of memory issues.
Gene& Genome::getGeneByIndex(unsigned index, bool simulated)
{
	if (simulated)
	{
		bool checker = checkIndex(index, 1, (unsigned)simulatedGenes.size());

		if (!checker)
		{
			my_printError("Warning: Invalid index given for simulated genes, returning simulated gene 1.\n");
		}

		return checker ? simulatedGenes[index - 1] : simulatedGenes[0];
	}
	else
	{
		bool checker = checkIndex(index, 1, (unsigned)genes.size());

		if (!checker)
		{
			my_printError("Warning: Invalid index given for genes, returning gene 1.\n");
		}

		return checker ? genes[index - 1] : genes[0];
	}
}


Gene& Genome::getGeneById(std::string ID, bool simulated)
{
	return getGene(ID, simulated);
}


Genome Genome::getGenomeForGeneIndicesR(std::vector <unsigned> indices, bool simulated)
{
	Genome genome;

	for (unsigned i = 0; i < indices.size(); i++)
	{
		if (indices[i] < 1 || indices[i] > getGenomeSize(simulated))
		{
			my_printError("Error in Genome::getGenomeForGeneIndices. An index specified is out of bounds for the genome!");
			my_printError("Returning empty Genome.");
			genome.clear();
			return genome;
		}
		else
		{
			simulated ? genome.addGene(simulatedGenes[indices[i]], true) : genome.addGene(genes[indices[i]], false);
		}
	}

	return genome;
}





//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//


RCPP_EXPOSED_CLASS(Gene)
RCPP_EXPOSED_CLASS(Genome)

RCPP_MODULE(Genome_mod)
{
	class_<Genome>("Genome")
		//Constructors & Destructors:
		.constructor("empty constructor")


		//File I/O Functions:
		.method("readFasta", &Genome::readFasta, "reads a genome into the object")
		.method("writeFasta", &Genome::writeFasta, "writes the genome to a fasta file")
		.method("readRFPFile", &Genome::readRFPFile, "reads RFP data in for the RFP model")
		.method("writeRFPFile", &Genome::writeRFPFile)
		.method("readObservedPhiValues", &Genome::readObservedPhiValues)


		//Gene Functions:
		.method("addGene", &Genome::addGene) //TEST THAT ONLY!
		.method("getGenes", &Genome::getGenes) //TEST THAT ONLY!
		.method("getNumGenesWithPhi", &Genome::getNumGenesWithPhi) //TEST THAT ONLY!


		//Other Functions:
		.method("checkIndex", &Genome::checkIndex) //TEST THAT ONLY!
		.method("getGenomeSize", &Genome::getGenomeSize, "returns how many genes are in the genome")
		.method("clear", &Genome::clear, "clears the genome")
		.method("getCodonCountsPerGene", &Genome::getCodonCountsPerGene,
			"returns a vector of codon counts for a given gene")



		//R Section:
		.method("getGeneByIndex", &Genome::getGeneByIndex, "returns a gene for a given index")
		.method("getGeneById", &Genome::getGeneById) //TEST THAT ONLY!
		.method("getGenomeForGeneIndices", &Genome::getGenomeForGeneIndicesR,
			"returns a new genome based on the ones requested in the given vector")
		;
}
#endif
