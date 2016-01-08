
#ifndef STANDALONE
#include "include/base/Trace.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )
  
    //Getter Functions:
    .method("getStdDevSynthesisRateAcceptanceRatioTrace", &Trace::getStdDevSynthesisRateAcceptanceRatioTrace)
    .method("getSynthesisRateTrace", &Trace::getSynthesisRateTrace)
    .method("getSynthesisRateAcceptanceRatioTrace", &Trace::getSynthesisRateAcceptanceRatioTrace)
    .method("getCodonSpecficAcceptanceRatioTraceForAA", &Trace::getCodonSpecficAcceptanceRatioTraceForAA)
    .method("getMixutreAssignmentTrace", &Trace::getMixtureAssignmentTrace)
    .method("getCodonSpecificAcceptanceRatioTrace", &Trace::getCodonSpecificAcceptanceRatioTrace)
    .method("getMixtureProbabilitiesTrace", &Trace::getMixtureProbabilitiesTrace)
    .method("getExpectedSynthesisRateTrace", &Trace::getExpectedSynthesisRateTrace)
    .method("getSynthesisOffsetAcceptanceRatioTrace", &Trace::getSynthesisOffsetAcceptanceRatioTrace)
    .method("getSynthesisOffsetAcceptanceRatioTraceForIndex", &Trace::getSynthesisOffsetAcceptanceRatioTraceForIndex)
    .method("getCodonSpecificParameterTrace", &Trace::getCodonSpecificParameterTrace)
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
    

    //ROC Specific:
    .method("getCodonSpecificParameterTraceByMixtureElementForCodon", &Trace::getCodonSpecificParameterTraceByMixtureElementForCodonR)
    .method("getSynthesisOffsetTrace", &Trace::getSynthesisOffsetTraceR)
    .method("getObservedSynthesisNoiseTrace", &Trace::getObservedSynthesisNoiseTraceR)
    .method("setSynthesisOffsetTrace", &Trace::setSynthesisOffsetTrace)
    .method("setSynthesisOffsetAcceptanceRatioTrace", &Trace::setSynthesisOffsetAcceptanceRatioTrace)
    .method("setObservedSynthesisNoiseTrace", &Trace::setObservedSynthesisNoiseTrace)
    .method("setCodonSpecificParameterTrace", &Trace::setCodonSpecificParameterTrace)

    ;
}
#endif

