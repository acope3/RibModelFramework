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

	if(!(this->genes == other.genes)) { match = false; std::cerr<<"fail1\n";} //Do a ! operation because only the gene comparison is implemented,
	//not the != operator.
	if(!(this->simulatedGenes == other.simulatedGenes)) { match = false; std::cerr<<"fail2\n"; }
	if(this->numGenesWithPhi != other.numGenesWithPhi) { match = false; std::cerr<<"fail3\n";}

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
		{
			clear();
		}
		std::ifstream Fin;
		Fin.open(filename.c_str());
		if (Fin.fail())
		{
#ifndef STANDALONE
			Rf_error("Error in Genome::readFasta: Can not open Fasta file %s\n", filename.c_str());
#else
			std::cerr << "Error in Genome::readFasta: Can not open Fasta file " << filename << "\n";
#endif
		}
		else
		{

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
					 2. Sequence line, which starts with any charactor except for ">"
					 3. End of file, detected by function eof().
				 */

				if(buf[0] == '>' ) newLine=1;
				else if(Fin.eof()) newLine=3;
				else newLine=2;

				if( newLine == 1 )
				{ // this is a start of a new chain.
					if(!fastaFormat)
					{
						// if it is the first chain, just build a new chain object
						tmpGene.clear();
						fastaFormat = true;
					} else
					{
						// otherwise, need to store the old chain first, then build a new chain
						tmpGene.setSequence(tempSeq);
						//tmpGene.cleanSeq();
						addGene(tmpGene);

						tmpGene.clear();
						tempSeq = "";
					}
					tmpGene.setDescription( buf.substr(1,buf.size()-1) );
					std::size_t pos = buf.find(" ") - 1;
					tmpGene.setId( buf.substr(1,pos) );
				}

				if( newLine == 2 && fastaFormat )
				{ // sequence line
					tempSeq.append(buf);
					//tmpGene.seq.append(buf);
				}

				if( newLine == 3 )
				{ // end of file
					if( !fastaFormat ) throw std::string("Genome::readFasta throws: ") + std::string(filename) + std::string(" is not in Fasta format.");
					else
					{
						// otherwise, need to store the old chain first, then to
						// build a new chain
						tmpGene.setSequence(tempSeq);
						//tmpGene.cleanSeq();
						addGene(tmpGene);
						break;
					}
				}
			} // end while
		} // end else
	} // end try
	catch(char* pMsg)
	{
#ifndef STANDALONE
		Rf_error("\nException: %s\n", pMsg);
#else
		std::cerr << std::endl << "Exception:" << pMsg << std::endl;
#endif
	}
}


void Genome::writeFasta (std::string filename, bool simulated)
{
	try {
		std::ofstream Fout;
		Fout.open(filename.c_str());
		if (Fout.fail())
		{
#ifndef STANDALONE
			Rf_error("Error in Genome::writeFasta: Can not open output Fasta file %s\n", filename.c_str());
#else
			std::cerr << "Error in Genome::writeFasta: Can not open output Fasta file " << filename << "\n";
#endif
		}
		else
		{
			if (simulated)
			{
				for(unsigned i = 0u; i < simulatedGenes.size(); i++)
				{
					Fout << ">" << simulatedGenes[i].getDescription() <<"\n";
					for(unsigned j = 0u; j < simulatedGenes[i].length(); j++)
					{
						Fout << simulatedGenes[i].getNucleotideAt(j);
						if((j + 1) % 60 == 0) Fout << std::endl;
					}
					Fout << std::endl;
				}
			}
			else
			{
				for(unsigned i = 0u; i < genes.size(); i++)
				{
					Fout << ">" << genes[i].getDescription() << std::endl;
					for(unsigned j = 0u; j < genes[i].length(); j++)
					{
						Fout << genes[i].getNucleotideAt(j);
						if((j + 1) % 60 == 0) Fout << std::endl;
					}
					Fout << std::endl;
				}
			}
		} // end else
		Fout.close();
	} // end try
	catch(char* pMsg)
	{
#ifndef STANDALONE
		Rf_error("\nException: %s\n", pMsg);
#else
		std::cerr << std::endl << "Exception:" << pMsg << std::endl;
#endif
	}
}

void Genome::readPANSEFile(std::string filename)
{
	std::ifstream Fin;
	Fin.open(filename.c_str());
	if (Fin.fail())
	{
#ifndef STANDALONE
		Rf_error("Error in Genome::readPANSEFile: Can not open PANSE file %s\n", filename.c_str());
#else
		std::cerr << "Error in Genome::readPANSEFile: Can not open PANSE file " << filename << "\n";
#endif
	}

	std::string tmp;
	std::getline(Fin, tmp); //trash the first line
	std::string prevID = "";
	Gene tmpGene;
	bool first = true;
	std::string seq = "";

	std::map <std::string, unsigned> genes;
	std::map <std::string, unsigned>::iterator git;

	//First pass: For each gene name, count its size.
	while (std::getline(Fin,tmp))
	{
		std::size_t pos = tmp.find(",");
		std::string ID = tmp.substr(0, pos);

		if (first)
		{
			prevID = ID;
			first = false;
			genes[ID] = 0;
		}
		if (ID != prevID)
		{
			genes[prevID] *= 3; //multiply by three since codons
			//std::cerr << "size of gene " << prevID << " is " << genes[prevID] << "\n";
			genes[ID] = 0; //initialize the new genes[ID]
		}
		prevID = ID;
		genes[ID]++;
	}
	genes[prevID] *= 3; //multiply by three since codons
	//std::cerr << "size of gene " << prevID << " is " << genes[prevID] << "\n";
	Fin.close();

	//std::cerr << "Re-traversing\n";

	Fin.open(filename.c_str());
	std::getline(Fin, tmp); //retrash first line
	prevID = "";
	first = true;

	//Now for each line associated with a gene ID, set the string appropriately
	while (std::getline(Fin,tmp))
	{
		std::size_t pos = tmp.find(",");
		std::string ID = tmp.substr(0, pos);

		if (first)
		{
			prevID = ID;
			first = false;
			git = genes.find(ID);
			seq.resize(git->second);
		}
		if (ID != prevID)
		{
			tmpGene.setId(prevID);
			tmpGene.setDescription("No description for PANSE Model");
			//tmpGene.setDescription(prevID + " Test Gene");
			tmpGene.setSequence(seq);
			addGene(tmpGene); //add to genome
			tmpGene.clear();
			std::cerr << prevID << " " << seq << "\n";
			seq = "";

			git = genes.find(ID);
			seq.resize(git->second);
		}
		// PANSE file format: GeneID,Codon,Position (1-indexed),rfp_count

		std::size_t pos2 = tmp.find(",", pos + 1);

		// Each codon is guaranteed to be of size 3
		std::string codon = tmp.substr(pos + 1, 3);

		// Make pos find next comma to find position integer
		pos = tmp.find(",", pos2 + 1);
		unsigned position = (unsigned)std::atoi(tmp.substr(pos2 + 1, pos - (pos2 + 1)).c_str()) - 1;
		//std::cerr << "position is " << position << " for codon " << codon << "\n";

		// Get rfp_count from last comma-separated value
		unsigned rfp_count = (unsigned)std::atoi(tmp.substr(pos + 1).c_str());

		// Now add codon to appropriate place in sequence
		seq[position * 3] = codon[0];
		seq[position * 3 + 1] = codon[1];
		seq[position * 3 + 2] = codon[2];

		prevID = ID;

		// Not sure what these do in context of PANSE
		//unsigned index = SequenceSummary::codonToIndex(codon);
		//tmpGene.geneData.setRFPObserved(index, tmpRFP);
	}

	tmpGene.setId(prevID);
	tmpGene.setDescription("No description for PANSE Model");
	//tmpGene.setDescription(prevID + " Test Gene");
	tmpGene.setSequence(seq);

	addGene(tmpGene); //add to genome
	std::cerr << prevID << " " << seq << "\n";

	Fin.close();
}


void Genome::readRFPFile(std::string filename)
{
	std::ifstream Fin;
	Fin.open(filename.c_str());
	if (Fin.fail())
	{
#ifndef STANDALONE
			Rf_error("Error in Genome::readRFPFile: Can not open RFP file %s\n", filename.c_str());
#else
			std::cerr << "Error in Genome::readRFPFile: Can not open RFP file " << filename << "\n";
#endif
	}

	std::string tmp;
	std::getline(Fin, tmp); //trash the first line
	std::string prevID = "";
	Gene tmpGene;
	bool first = true;
	std::string seq = "";

	while (std::getline(Fin,tmp))
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
			//tmpGene.setDescription("No description for RFP Model");
			tmpGene.setSequence(seq);
			addGene(tmpGene); //add to genome
			tmpGene.clear();
			seq = "";
		}
		std::size_t pos2 = tmp.find(",", pos + 1);
		std::string value = tmp.substr(pos + 1, pos2 - (pos + 1));
		unsigned tmpRFP = (unsigned)std::atoi(value.c_str());
		pos = tmp.find(",", pos2 + 1);
		value = tmp.substr(pos2 + 1, pos - (pos2 + 1));
		unsigned counts = (unsigned)std::atoi(value.c_str());

		std::string codon = tmp.substr(pos + 1, 3);
		for (unsigned i = 0; i < counts; i++)
			seq += codon;

		prevID = ID;
		unsigned index = SequenceSummary::codonToIndex(codon);
		tmpGene.geneData.setRFPObserved(index, tmpRFP);
	}

	tmpGene.setId(prevID);
	tmpGene.setDescription("No description for RFP Model");
	tmpGene.setSequence(seq);

	addGene(tmpGene); //add to genome

	Fin.close();
}


void Genome::writeRFPFile(std::string filename, bool simulated)
{
	std::ofstream Fout;
	Fout.open(filename.c_str());
	if (Fout.fail())
	{
#ifndef STANDALONE
		Rf_error("Error in Genome::writeRFPFile: Can not open output RFP file %s\n", filename.c_str());
#else
		std::cerr <<"Error in Genome::writeRFPFile: Can not open output RFP file " << filename <<"\n";
#endif
	}

	Fout <<"ORF,RFP_Counts,Codon_Counts,Codon\n";
	for (unsigned geneIndex = 0; geneIndex < genes.size(); geneIndex++)
	{
		Gene *currentGene;
		if (simulated)
			currentGene = &simulatedGenes[geneIndex];
		else
			currentGene = &genes[geneIndex];


		for (unsigned codonIndex = 0; codonIndex < 64; codonIndex++)
		{
			std::string codon = SequenceSummary::codonArray[codonIndex];

			Fout << currentGene->getId() <<",";
			Fout << currentGene->geneData.getRFPObserved(codonIndex) <<",";
			Fout << currentGene->geneData.getCodonCountForCodon(codonIndex) <<"," << codon <<"\n";
		}
	}
	Fout.close();
}


void Genome::readObservedPhiValues(std::string filename, bool byId)
{
	std::ifstream input;
	std::string tmp;
	unsigned numPhi = 0;
	bool exitfunction = false;

	input.open(filename);
	if (input.fail())
	{
#ifndef STANDALONE
		Rf_error("Error in Genome::readObservedPhiValues: Can not open file %s\n", filename.c_str());
#else
		std::cerr << "Error in Genome::readObservedPhiValues: Can not open file " << filename << "\n";
#endif
	} //End of opening a file and it resulting in a failure.
	else
	{
		std::getline(input, tmp); //Trash the header line

		if (genes.size() == 0)
		{
#ifndef STANDALONE
			Rf_error("Genome is empty, function will not execute!\n");
#else
			std::cerr << "Genome is empty, function will not execute!\n";
#endif
		} //end of checking size constraint.
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
					{
#ifndef STANDALONE
                        Rf_warning("Gene %d not found!\n", geneID.c_str());
#else
                        std::cerr << "Gene " << geneID << " not found!\n";
#endif
                    }
                    else //gene is found
                    {
                        unsigned phiGrouping = 0;
                        std::string val = "";
                        bool notDone = true;
						double dval;

                        while (notDone)
						{
							std::size_t pos2 = tmp.find(",", pos + 1);

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
#ifndef STANDALONE
								Rf_warning("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.");
#else
								std::cerr <<
								"WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.\"";
#endif
							}
							it->second->observedSynthesisRateValues.push_back(dval);
						}

						if (first)
						{
							first = false;
							numPhi = (unsigned) it->second->observedSynthesisRateValues.size();
							numGenesWithPhi.resize(numPhi, 0);
						}
						else if (it->second->observedSynthesisRateValues.size() != numGenesWithPhi.size())
						{
#ifndef STANDALONE
							Rf_error("Gene %d: has a different number of phi values given other genes: \n", geneID.c_str());
                            Rf_error("%d\n", it-> second -> getId().c_str());
                            Rf_error("%d\n", it-> second -> observedSynthesisRateValues.size() + 1);
                            Rf_error("Exiting function.\n");
#else
							std::cerr << geneID << ": has a different number of phi values given than other genes: ";
							std::cerr << it->second->getId() << "\n";
							std::cerr << it->second->observedSynthesisRateValues.size() + 1 << ". Exiting function.\n";
#endif
							exitfunction = true;
							for (unsigned a = 0; a < getGenomeSize(); a++) {
								genes[a].observedSynthesisRateValues.clear();
							}
							break;
						}
                    }
                    if (exitfunction) break;
                }
                if (!exitfunction)
				{
                    for (unsigned geneIndex = 0; geneIndex < getGenomeSize(); geneIndex++)
					{
						Gene *gene = &(getGene(geneIndex));
                        if (getGene(geneIndex).observedSynthesisRateValues.size() != numPhi)
						{
#ifndef STANDALONE
                            Rf_warning("Gene # %d (%s) does not have any phi values.", geneIndex, gene->getId().c_str());
                            Rf_warning(" Filling with -1's\n");
#else
                            std::cerr << "Gene # " << geneIndex << " (" << gene->getId() <<
                            ") does not have any phi values.";
                            std::cerr << " Filling with -1's\n";
#endif
                            gene->observedSynthesisRateValues.resize(numPhi, -1);
                        }

						for (unsigned i = 0; i < numPhi; i++)
						{
							if (gene->observedSynthesisRateValues[i] != -1)
							{
								numGenesWithPhi[i]++;
								std::cerr << "For index i = " << i << " numGenesWithPhi[i] is " << numGenesWithPhi[i]
									<< "\n";
							}
						}
                    }
                }
            }
                //
                //
                //end of putting in by ID
			    //
                //
            else //doing this by index
			{
				unsigned geneIndex = 0;
				bool first = true;
				bool firstValue = false;
				while (std::getline(input, tmp))
				{
					// HOLLIS TODO: Let's change it to a vector-based setup like
					// In byID.

					//std::size_t pos = tmp.find(",");


					// Begin old code

					if (geneIndex >= genes.size())
					{
#ifndef STANDALONE
						Rf_error("GeneIndex exceeds the number of genes in the genome. Exiting function\n");
#else
						std::cerr << "GeneIndex exceeds the number of genes in the genome. Exiting function\n";
#endif
						break;
					}
					std::size_t pos = tmp.find(",");

					std::string val = "";
					bool notDone = true;
					double dval;

					while (notDone)
					{
						std::size_t pos2 = tmp.find(",", pos + 1);

						if (pos2 == std::string::npos)
						{
							val = tmp.substr(pos + 1, tmp.size() + 1);
							notDone = false;

							//I think this code is in the wrong place. Moving down.

							if (first)
							{
								first = false;
								numPhi = (unsigned) genes[geneIndex].observedSynthesisRateValues.size() + 1;
								numGenesWithPhi.resize(numPhi);
								//What does this do? It just means first value
								//Is initialized at 1 more than the rest.
								//if (firstValue) numGenesWithPhi[0] = 1;
							}
							else if (numPhi != genes[geneIndex].observedSynthesisRateValues.size() + 1)
							{
#ifndef STANDALONE
                                Rf_error("Gene %d: has a different number of phi values given other genes: \n", geneIndex);
                                Rf_error("%d\n", genes[geneIndex].observedSynthesisRateValues.size() + 1);
                                Rf_error("Exiting function.\n");
#else
								std::cerr << "Gene " << geneIndex <<
								": has a different number of phi values given than other genes: \n";
								std::cerr << genes[geneIndex].observedSynthesisRateValues.size() + 1 <<
								". Exiting function\n";
#endif
								exitfunction = true;
								for (unsigned a = 0; a < getGenomeSize(); a++)
								{
									genes[a].observedSynthesisRateValues.clear();
								}
								break;
							}

						}
						else
						{
							val = tmp.substr(pos + 1, pos2 - (pos + 1));
							//pos = pos2;
						}
						dval = std::atof(val.c_str());

						//If the value is negative, nan, or 0, we set it to -1 and throw a warning message.
						if (dval <= 0 || std::isnan(dval))
						{
							dval = -1;
#ifndef STANDALONE
							Rf_warning("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.");
#else
							std::cerr <<
							"WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Negative Value stored.\"";
#endif
						}
						//genes[geneIndex].observedSynthesisRateValues.push_back(dval);
					//}

					//HOLLIS TODO: is this shit even right
					/*
					if (first) {
						first = false;
						numPhi = (unsigned) genes[geneIndex].observedSynthesisRateValues.size() + 1;
						numGenesWithPhi.resize(numPhi);
						//What does this do? It just means first value
						//Is initialized at 1 more than the rest.
						//if (firstValue) numGenesWithPhi[0] = 1;
					}
					else if (numPhi != genes[geneIndex].observedSynthesisRateValues.size() + 1) {
#ifndef STANDALONE
						Rf_error("Gene %d: has a different number of phi values given other genes: \n", geneIndex);
                                Rf_error("%d\n", genes[geneIndex].observedSynthesisRateValues.size() + 1);
                                Rf_error("Exiting function.\n");
#else
						std::cerr << "Gene " << geneIndex <<
						": has a different number of phi values given than other genes: \n";
						std::cerr << genes[geneIndex].observedSynthesisRateValues.size() + 1 <<
						". Exiting function\n";
#endif
						exitfunction = true;
						for (unsigned a = 0; a < getGenomeSize(); a++) {
							genes[a].observedSynthesisRateValues.clear();
						}
						break;
					} else{
						numGenesWithPhi[genes[geneIndex].observedSynthesisRateValues.size()]++;
					}
					 */

						else
						{
							if (first)
							{
								firstValue = true;
							}
							else
							{
								numGenesWithPhi[genes[geneIndex].observedSynthesisRateValues.size()]++;
							}
						}

						genes[geneIndex].observedSynthesisRateValues.push_back(dval);
						pos = pos2;
					}

					geneIndex++;
					if (exitfunction) break;
				}

				//POSSIBLE TODO: Parallelize this with byID's code for readability.
				//Current syntax is ambiguous.
				if (!exitfunction)
				{
					for (unsigned i = 0; i < numPhi; i++)
					{

						if (numGenesWithPhi[i] != getGenomeSize())
						{
#ifndef STANDALONE
                            Rf_warning("The last %d genes do not have phi values.", getGenomeSize() - numGenesWithPhi[i]);
                            Rf_warning("Please check your file to make sure every gene has a phi value. Filling empty genes ");
                            Rf_warning("with -1's for calculations.\n");
#else
							std::cerr << "The last " << getGenomeSize() - numGenesWithPhi[i] <<
							" genes do not have phi values.";
							std::cerr <<
							"Please check your file to make sure every gene has a phi value. Filling empty genes ";
							std::cerr << "with -1's for calculations.\n";
#endif

							for (unsigned a = numGenesWithPhi[i]; a < getGenomeSize(); a++)
							{
								Gene *gene = &(getGene(a));
								gene->observedSynthesisRateValues[i] = -1;
							}
						}

						//Did I add this here?
						std::cerr << "For index i = " << i << " numGenesWithPhi[i] is " << numGenesWithPhi[i]
						<< "\n";
						//
					}
				}
			}//end of reading by index
		}
		input.close();
	}
}





//------------------------------------//
//---------- Gene Functions ----------//
//------------------------------------//


void Genome::addGene(const Gene& gene, bool simulated)
{
	if (!simulated)
		genes.push_back(gene);
	else
		simulatedGenes.push_back(gene);
}


std::vector <Gene> Genome::getGenes(bool simulated)
{
	return  !simulated ? genes : simulatedGenes;
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
		if (!simulated) tempGene = genes[geneIndex];
		else tempGene = simulatedGenes[geneIndex];
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


Genome Genome::getGenomeForGeneIndicies(std::vector <unsigned> indicies, bool simulated)
{
	Genome genome;

	for (unsigned i = 0; i < indicies.size(); i++)
	{
		simulated ? genome.addGene(simulatedGenes[indicies[i]]) : genome.addGene(genes[indicies[i]]);
	}

	return genome;
}


std::vector<unsigned> Genome::getCodonCountsPerGene(std::string codon)
{
	std::vector<unsigned> codonCounts(genes.size());
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	for(unsigned i = 0u; i < genes.size(); i++)
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


std::vector<unsigned> Genome::getNumGenesWithPhi()
{
	return numGenesWithPhi;
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
		Rf_error("Index: %d is out of bounds. Index must be between %d & %d\n", index, lowerbound, upperbound);
	}
	return check;
}


Gene& Genome::getGeneByIndex(unsigned index, bool simulated) //NOTE: This function does the check and performs the function itself because of memory issues.
{
	bool checker = checkIndex(index, 1, (unsigned)genes.size());
	if (!checker)
	{
#ifndef STANDALONE
		Rf_warning("Invalid index given, returning gene 1, not simulated\n");
#else
		std::cerr << "Invalid index given, returning gene 1, not simulated\n";
#endif
	}
	return checker ? simulated ? simulatedGenes[index - 1] : genes[index - 1] : genes[0];
}


Gene& Genome::getGeneById(std::string ID, bool simulated)
{
	return getGene(ID, simulated);
}


Genome Genome::getGenomeForGeneIndiciesR(std::vector <unsigned> indicies, bool simulated)
{
	Genome genome;
	bool check = true;
	for (unsigned i = 0; i < indicies.size(); i++)
	{
		if (indicies[i] < 1 || indicies[i] > getGenomeSize())
		{
			check = false;
			break;
		}
		else
		{
			indicies[i] -= 1;
		}
	}
	return check ? getGenomeForGeneIndicies(indicies, simulated) : genome;
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
		.method("getCodonCountsPerGene", &Genome::getCodonCountsPerGene, "returns a vector of codon counts for a given gene")



		//R Section:
		.method("getGeneByIndex", &Genome::getGeneByIndex, "returns a gene for a given index")
		.method("getGeneById", &Genome::getGeneById) //TEST THAT ONLY!
		.method("getGenomeForGeneIndicies", &Genome::getGenomeForGeneIndiciesR, "returns a new genome based on the ones requested in the given vector")
		;
}
#endif
