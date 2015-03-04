#include "../include/ROCModel.h"

#include <vector>
#include <math.h>

ROCModel::ROCModel()
{
    //ctor
}

ROCModel::~ROCModel()
{
    //dtor
}

ROCModel::ROCModel(const ROCModel& other)
{
    //copy ctor
}

double ROCModel::calculateLogLiklihoodPerGene(Gene& gene, int geneIndex, ROCParameter& parameter, bool proposed)
{
    double logLikelihood = 0.0;
    double phiValue = parameter.getExpression(geneIndex, proposed);

    SequenceSummary seqsum = gene.getSequenceSummary();
    for(int i = 0; i < 22; i++)
    {
        char curAA = seqsum.AminoAcidArray[i];
        // skip amino acids with only one codon or stop codons
        if(curAA == 'X' || curAA == 'M' || curAA == 'W') continue;
        // skip amino acids which do not occur in current gene. Avoid useless calculations and multiplying by 0
        if(seqsum.getAAcountForAA(i) == 0) continue;

        // get codon count (total count not parameter count)
        int numCodons = seqsum.GetNumCodonsForAA(curAA);
        int codonCount[numCodons];
        // prepare array for codon counts for AA
        unsigned* codonRange = SequenceSummary::AAToCodonRange(curAA);
        // get codon counts for AA
        unsigned j = 0u;
        for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
        {
            codonCount[j] = seqsum.getCodonCountForCodon(i);
        }
        // get mutation and selection parameter for gene
        double mutation[numCodons - 1];
        parameter.getParameterForCategory(gene.getMutationCategory(), ROCParameter::dM, curAA, false, mutation);
        double selection[numCodons - 1];
        parameter.getParameterForCategory(gene.getDeltaEtaCategory(), ROCParameter::dEta, curAA, false, selection);

        // calculate codon probabilities
        double codonProbabilities[numCodons];
        calculateCodonProbabilityVector(numCodons, mutation, selection, phiValue, codonProbabilities);

        // calculate likelihood for current AA for this combination of selection and mutation category
        for(unsigned i = 0; i < numCodons; i++)
        {
            logLikelihood += std::log(codonProbabilities[i]) * codonCount[i];
        }
    }
    double sPhi = parameter.getSphi(false);
    return logLikelihood - std::log(ROCParameter::densityLogNorm(phiValue, -(sPhi * sPhi) / 2, sPhi));
}

void ROCModel::calculateCodonProbabilityVector(unsigned numCodons, double mutation[], double selection[], double phi, double codonProb[])
{
    // calculate numerator and denominator for codon probabilities
    codonProb[0] = 1.0;
    double denominator = 1.0;
    for(unsigned i = 1; i < numCodons; i++)
    {
        codonProb[i] = std::exp( -mutation[i - 1] - (selection[i - 1] * phi) );
        denominator += codonProb[i];
    }
    // normalize codon probabilities
    for(unsigned i = 0; i < numCodons; i++) { codonProb[i] = codonProb[i] / denominator;}
}

/*
double ROCModel::calculateLogLikelihood(Genome genome, ROCParameter parameter, unsigned category, unsigned paramType, bool proposed)
{
    int ngenes = genome.getGenomeSize();
    double likelihood = 0.0;
    for(int i = 0; i < ngenes; i++)
    {
        likelihood += calculateLikelihoodPerGene(genome.getGene(i), parameter, proposed);
    }
    return likelihood;
}



double ROCModel::calculateLikelihoodPerAAForCategory(char aa, Genome& genome, ROCParameter& parameter, unsigned paramType, unsigned category)
{

}

double ROCModel::calculateLikelihoodPerAA(char aa, Genome& genome, ROCParameter& parameter, bool pdM, bool pdEta)
{
    double logLikelihood = 0.0;
    unsigned numCodons = SequenceSummary::GetNumCodonsForAA(aa);

    // prepare array for codon counts for AA
    unsigned* codonRange = SequenceSummary::AAToCodonRange(aa);
    int codonCount[numCodons];

    // iterate over all genes
    unsigned nGenes = genome.getGenomeSize();
    for(unsigned curGeneIdx = 0; curGeneIdx < nGenes; curGeneIdx++)
    {
        double geneLogLikelihood = 0.0;
        Gene gene = genome.getGene(curGeneIdx);
        SequenceSummary seqsum = gene.getSequenceSummary();

        // get codon counts for AA
        unsigned j = 0u;
        for(unsigned i = codonRange[0]; i < codonRange[1]; i++, j++)
        {
            codonCount[j] = seqsum.getCodonCountForCodon(i);
        }

        // iterate over all categories
        std::vector<Category> etaCat = gene.getDeltaEtaCategories();
        // all selection categories
        for(unsigned etaCategoryCount = 0; etaCategoryCount < etaCat.size(); etaCategoryCount++)
        {
            double etaCategoryLogLik = 0.0;
            double* selection = parameter.getParameterForCategory(etaCat[etaCategoryCount].cat, ROCParameter::dEta, aa, pdEta);
            std::vector<Category> mCat = gene.getMutationCategories();
            // all mutation categories
            for(unsigned mCategoryCount = 0; mCategoryCount < mCat.size(); mCategoryCount++)
            {
                double mCategoryLogLik = 0.0;
                double* mutation = parameter.getParameterForCategory(mCat[mCategoryCount].cat, ROCParameter::dM, aa, pdM);
                // calculate codon probability vector
                double* codonProb = calculateCodonProbabilityVector(numCodons, mutation, selection, gene.getPhi());
                // calculate likelihood for current AA for this combination of selection and mutation category
                for(unsigned i = 0; i < numCodons; i++) { mCategoryLogLik += codonProb[i] * codonCount[i];}
                // sum over all mutation categories
                etaCategoryLogLik += mCat[mCategoryCount].prob * mCategoryLogLik;
            }
            // sum over all selection categories
            geneLogLikelihood += etaCat[etaCategoryCount].prob * etaCategoryLogLik;
        }
        // sum over all genes
        logLikelihood += geneLogLikelihood;
    }
    return logLikelihood;
}

double ROCModel::calculateLikelihoodPerGene(Gene gene, ROCParameter parameter, bool propsed)
{
    SequenceSummary seqsum = gene.getSequenceSummary();
    double logLikelihood = 0.0;
//    for(unsigned categoryCount = 0; categoryCount < gene.categories.size(); categoryCount++)
//    {
//        // calculate likelihood for gene i in category j
//        double categoryLogLik = 0.0;
//        unsigned startAAIndex = 0;
//        for(int aaIndex = 0; aaIndex < 22; aaIndex++)
//        {
//            char aa = SequenceSummary::AminoAcidArray[aaIndex];
//            if(aa == 'M' || aa == 'W') continue; //ignore one codon AA
//
//            // to find position of parameters for codon for AA
//            unsigned numCodons = seqsum.GetNumCodonsForAA(aa);
//            unsigned endAAIndex = numCodons + startAAIndex;
//
//            // get parameters for the category this gene is in.
//            std::vector<double> mutation = parameter.getParameterForCategory(gene.categories[categoryCount].cat, ROCParameter::dM);
//            std::vector<double> selection = parameter.getParameterForCategory(gene.categories[categoryCount].cat, ROCParameter::dEta);
//            // get genes phi values
//            double phi = gene.getPhi();
//
//            double logCodonProp[numCodons];
//            int codonCount[numCodons];
//            double denominator = 0.0;
//
//            // calculate numerator and denominator for codon probabilities
//            for(unsigned i = 0; startAAIndex < endAAIndex; startAAIndex++, i++)
//            {
//                codonCount[i] = seqsum.getCodonCountForCodon(startAAIndex);
//                logCodonProp[i] = mutation[startAAIndex] + selection[startAAIndex] * phi;
//                denominator += codonProp[i];
//            }
//
//            // The next two loops can be combined into one. I left it as two for code readability!
//            // normalize codon probabilities
//            for(unsigned i = 0; i < numCodons; i++) { logCodonProp[i] = logCodonProp[i] - denominator;}
//            // calculate likelihood for this codon in this category, sum over all codons for AA and all AA in gene.
//            for(unsigned i = 0; i < numCodons; i++) { categoryLogLik += std::exp(logCodonProp[i]) * codonCount[i];}
//        }
//        // Expected likelihood over all categories the gene is in.
//        logLikelihood += gene.categories[categoryCount].prob * categoryLogLik;
//    }
    return logLikelihood;
}

*/




