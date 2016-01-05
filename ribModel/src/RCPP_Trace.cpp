#ifndef STANDALONE
#include "include/base/Trace.h"
#include <Rcpp.h>
using namespace Rcpp;


RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )

    //Getter Functions:
    .method("getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene", &Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR)
    .method("getSynthesisRateTraceForGene", &Trace::getSynthesisRateTraceForGeneR)
    .method("getSynthesisRateTraceByMixtureElementForGene", &Trace::getSynthesisRateTraceByMixtureElementForGeneR)
    .method("getMixtureAssignmentTraceForGene", &Trace::getMixtureAssignmentTraceForGeneR)
    .method("getMixtureProbabilitiesTraceForMixture", &Trace::getMixtureProbabilitiesTraceForMixtureR)
    .method("getStdDevSynthesisRateTraces", &Trace::getStdDevSynthesisRateTraces)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)



    //Setter Functions:
    .method("setStdDevSynthesisRateTraces", &Trace::setStdDevSynthesisRateTraces)
    .method("setStdDevSynthesisRateAcceptanceRatioTrace", &Trace::setStdDevSynthesisRateAcceptanceRatioTrace)
    .method("setSynthesisRateTrace", &Trace::setSynthesisRateTrace)
    .method("setSynthesisRateAcceptanceRatioTrace", &Trace::setSynthesisRateAcceptanceRatioTrace)
    .method("setMixtureAssignmentTrace", &Trace::setMixtureAssignmentTrace)
    .method("setMixtureProbabilitiesTrace", &Trace::setMixtureProbabilitiesTrace)
    .method("setCodonSpecificAcceptanceRatioTrace", &Trace::setCodonSpecificAcceptanceRatioTrace)
    //.method("setCategories", &Trace::setCategories) //TODO: used?



    //ROC Specific:
	.method("getCodonSpecificParameterTraceByMixtureElementForCodon", &Trace::getCodonSpecificParameterTraceByMixtureElementForCodon)






    .method("getMixutreAssignmentTrace", &Trace::getMixtureAssignmentTrace)
    .method("getStdDevSynthesisRateAcceptanceRatioTrace", &Trace::getStdDevSynthesisRateAcceptanceRatioTrace)
    .method("getCodonSpecficAcceptanceRatioTraceForAA", &Trace::getCodonSpecficAcceptanceRatioTraceForAA)
    .method("getCodonSpecificAcceptanceRatioTrace", &Trace::getCodonSpecificAcceptanceRatioTrace)


    //These methods have specific R wrappers
    .method("getSynthesisRateAcceptanceRatioTrace", &Trace::getSynthesisRateAcceptanceRatioTrace)
    .method("getSynthesisRateTrace", &Trace::getSynthesisRateTrace)
    .method("getMixtureProbabilitiesTrace", &Trace::getMixtureProbabilitiesTrace)
    .method("getExpectedSynthesisRateTrace", &Trace::getExpectedSynthesisRateTrace)


	.method("getSynthesisOffsetAcceptanceRatioTrace", &Trace::getSynthesisOffsetAcceptanceRatioTrace)
	.method("getSynthesisOffsetTrace", &Trace::getSynthesisOffsetTrace)
	.method("getSynthesisOffsetAcceptanceRatioTraceForIndex", &Trace::getSynthesisOffsetAcceptanceRatioTraceForIndex)
	.method("getObservedSynthesisNoiseTrace", &Trace::getObservedSynthesisNoiseTrace)
	.method("setSynthesisOffsetTrace", &Trace::setSynthesisOffsetTrace)
	.method("setSynthesisOffsetAcceptanceRatioTrace", &Trace::setSynthesisOffsetAcceptanceRatioTrace)
	.method("setObservedSynthesisNoiseTrace", &Trace::setObservedSynthesisNoiseTrace)
	.method("setCodonSpecificParameterTrace", &Trace::setCodonSpecificParameterTrace)
	.method("getCodonSpecificParameterTrace", &Trace::getCodonSpecificParameterTrace)
    ;
}
#endif

