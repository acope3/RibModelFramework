
#ifndef STANDALONE
#include "include/base/Trace.h"
#include "include/ROC/ROCTrace.h"
#include "include/RFP/RFPTrace.h"
#include "include/FONSE/FONSETrace.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )
    //These methods have only a C++ implementation
    //.method("getSphiTrace", &Trace::getSphiTrace)
    .method("getSphiTraces", &Trace::getSphiTraces)
    //.method("getAPhiTrace", &Trace::getAPhiTrace)
    .method("getSphiAcceptanceRatioTrace", &Trace::getSphiAcceptanceRatioTrace)
    .method("getCspAcceptanceRatioTraceForAA", &Trace::getCspAcceptanceRatioTraceForAA)
    .method("getCspAcceptanceRatioTrace", &Trace::getCspAcceptanceRatioTrace)


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
    .method("getExpectedPhiTrace", &Trace::getExpectedPhiTrace)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)
    
    .method("setSphiTraces", &Trace::setSphiTraces)
    .method("setSphiAcceptanceRatioTrace", &Trace::setSphiAcceptanceRatioTrace)
    .method("setSynthesisRateTrace", &Trace::setSynthesisRateTrace)
    .method("setSynthesisRateAcceptanceRatioTrace", &Trace::setSynthesisRateAcceptanceRatioTrace)
    .method("setMixtureAssignmentTrace", &Trace::setMixtureAssignmentTrace)
    .method("setMixtureProbabilitiesTrace", &Trace::setMixtureProbabilitiesTrace)
    .method("setCspAcceptanceRatioTrace", &Trace::setCspAcceptanceRatioTrace)

    ;

  class_<ROCTrace>( "ROCTrace" )
    .derives<Trace>("Trace")
    //These methods have specific R wrappers
    .method("getMutationParameterTraceByMixtureElementForCodon", &ROCTrace::getMutationParameterTraceByMixtureElementForCodonR)
    .method("getSelectionParameterTraceByMixtureElementForCodon", &ROCTrace::getSelectionParameterTraceByMixtureElementForCodonR)
    .method("getMutationParameterTrace", &ROCTrace::getMutationParameterTrace)
    .method("getSelectionParameterTrace", &ROCTrace::getSelectionParameterTrace)
    .method("getAphiAcceptanceRatioTrace", &ROCTrace::getAphiAcceptanceRatioTrace)
	.method("getAphiTraces", &ROCTrace::getAphiTraceR)
	.method("getAphiAcceptanceRatioTraceForIndex", &ROCTrace::getAphiAcceptanceRatioTraceForIndex)
	.method("getSepsilonTraces", &ROCTrace::getSepsilonTraceR)
	.method("setAphiTrace", &ROCTrace::setAphiTrace)
    .method("setAphiAcceptanceRatioTrace", &ROCTrace::setAphiAcceptanceRatioTrace)
	.method("setSepsilonTrace", &ROCTrace::setSepsilonTrace)
	.method("setMutationParameterTrace", &ROCTrace::setMutationParameterTrace)
	.method("setSelectionParameterTrace", &ROCTrace::setSelectionParameterTrace)

    ;

	class_<RFPTrace>("RFPTrace")
		.derives<Trace>("Trace")
		.method("getAlphaParameterTraceByMixtureElementForCodon", &RFPTrace::getAlphaParameterTraceByMixtureElementForCodonR)
		.method("getAlphaParameterTrace", &RFPTrace::getAlphaParameterTrace)
		.method("getLambdaPrimeParameterTraceByMixtureElementForCodon", &RFPTrace::getLambdaPrimeParameterTraceByMixtureElementForCodonR)
		.method("getLambdaPrimeParameterTrace", &RFPTrace::getLambdaPrimeParameterTrace)
        .method("setAlphaParameterTrace", &RFPTrace::setAlphaParameterTrace)
        .method("setLambdaPrimeParameterTrace", &RFPTrace::setLambdaPrimeParameterTrace)
		;

	class_<FONSETrace>("FONSETrace")
		.derives<Trace>("Trace")

		.method("getMutationParameterTrace", &FONSETrace::getMutationParameterTrace)
		.method("getSelectionParameterTrace", &FONSETrace::getSelectionParameterTrace)
		//These methods have specific R wrappers
		.method("getMutationParameterTraceByMixtureElementForCodon", &FONSETrace::getMutationParameterTraceByMixtureElementForCodonR)
		.method("getSelectionParameterTraceByMixtureElementForCodon", &FONSETrace::getSelectionParameterTraceByMixtureElementForCodonR)
		;
}
#endif

