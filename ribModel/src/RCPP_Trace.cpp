
#ifndef STANDALONE
#include "include/base/Trace.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )
    //These methods have only a C++ implementation
    .method("getStdDevSynthesisRateTraces", &Trace::getStdDevSynthesisRateTraces)
    .method("getStdDevSynthesisRateAcceptanceRatioTrace", &Trace::getStdDevSynthesisRateAcceptanceRatioTrace)
    .method("getCodonSpecficAcceptanceRatioTraceForAA", &Trace::getCodonSpecficAcceptanceRatioTraceForAA)
    .method("getCodonSpecificAcceptanceRatioTrace", &Trace::getCodonSpecificAcceptanceRatioTrace)


    //These methods have specific R wrappers
    .method("getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene", &Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR)
    .method("getSynthesisRateAcceptanceRatioTrace", &Trace::getSynthesisRateAcceptanceRatioTrace)
    .method("getSynthesisRateTraceForGene", &Trace::getSynthesisRateTraceForGeneR)
    .method("getSynthesisRateTraceByMixtureElementForGene", &Trace::getSynthesisRateTraceByMixtureElementForGeneR)
    .method("getSynthesisRateTrace", &Trace::getSynthesisRateTrace)
    .method("getMixtureAssignmentTraceForGene", &Trace::getMixtureAssignmentTraceForGeneR)
    .method("getMixtureProbabilitiesTraceForMixture", &Trace::getMixtureProbabilitiesTraceForMixtureR)
    .method("getMixutreAssignmentTrace", &Trace::getMixtureAssignmentTrace)
    .method("getMixtureProbabilitiesTrace", &Trace::getMixtureProbabilitiesTrace)
    .method("getExpectedSynthesisRateTrace", &Trace::getExpectedSynthesisRateTrace)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)
    
    .method("setStdDevSynthesisRateTraces", &Trace::setStdDevSynthesisRateTraces)
    .method("setStdDevSynthesisRateAcceptanceRatioTrace", &Trace::setStdDevSynthesisRateAcceptanceRatioTrace)
    .method("setSynthesisRateTrace", &Trace::setSynthesisRateTrace)
    .method("setSynthesisRateAcceptanceRatioTrace", &Trace::setSynthesisRateAcceptanceRatioTrace)
    .method("setMixtureAssignmentTrace", &Trace::setMixtureAssignmentTrace)
    .method("setMixtureProbabilitiesTrace", &Trace::setMixtureProbabilitiesTrace)
    .method("setCodonSpecificAcceptanceRatioTrace", &Trace::setCodonSpecificAcceptanceRatioTrace)

	.method("getSynthesisOffsetAcceptanceRatioTrace", &ROCTrace::getSynthesisOffsetAcceptanceRatioTrace)
	.method("getSynthesisOffsetTrace", &ROCTrace::getAphiTraceR)
	.method("getSynthesisOffsetAcceptanceRatioTraceForIndex", &ROCTrace::getSynthesisOffsetAcceptanceRatioTraceForIndex)
	.method("getObservedSynthesisNoiseTrace", &ROCTrace::getObservedSynthesisNoiseTrace)
	.method("setSynthesisOffsetTrace", &ROCTrace::setSynthesisOffsetTrace)
	.method("setSynthesisOffsetAcceptanceRatioTrace", &ROCTrace::setSynthesisOffsetAcceptanceRatioTrace)
	.method("setObservedSynthesisNoiseTrace", &ROCTrace::setObservedSynthesisNoiseTrace)
	.method("setCodonSpecificParameterTrace", &ROCTrace::setCodonSpecificParameterTrace)
	.method("getCodonSpecificParameterTrace", &ROCTrace::getCodonSpecificParameterTrace)
	.method("getCodonSpecificParameterTraceByMixtureElementForCodon", &ROCTrace::getCodonSpecificParameterTraceByMixtureElementForCodon)
    ;
}
#endif

