
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


    //These methods have specific R wrappers
    .method("getSynthesisRateAcceptanceRatioTraceByMixtureElementForGene", &Trace::getSynthesisRateAcceptanceRatioTraceByMixtureElementForGeneR)
    .method("getSynthesisRateTraceForGene", &Trace::getSynthesisRateTraceForGeneR)
    .method("getSynthesisRateTraceByMixtureElementForGene", &Trace::getSynthesisRateTraceByMixtureElementForGeneR)
    .method("getMixtureAssignmentTraceForGene", &Trace::getMixtureAssignmentTraceForGeneR)
    .method("getMixtureProbabilitiesTraceForMixture", &Trace::getMixtureProbabilitiesTraceForMixtureR)
    .method("getExpectedPhiTrace", &Trace::getExpectedPhiTrace)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)
    ;

  class_<ROCTrace>( "ROCTrace" )
    .derives<Trace>("Trace")
    //These methods have specific R wrappers
    .method("getMutationParameterTraceByMixtureElementForCodon", &ROCTrace::getMutationParameterTraceByMixtureElementForCodonR)
    .method("getSelectionParameterTraceByMixtureElementForCodon", &ROCTrace::getSelectionParameterTraceByMixtureElementForCodonR)
	.method("getAphiTraces", &ROCTrace::getAphiTraceR)
	.method("getAphiAcceptanceRatioTrace", &ROCTrace::getAphiAcceptanceRatioTrace)
	.method("getSepsilonTraces", &ROCTrace::getSepsilonTraceR)
    ;

	class_<RFPTrace>("RFPTrace")
		.derives<Trace>("Trace")
		.method("getAlphaParameterTraceByMixtureElementForCodon", &RFPTrace::getAlphaParameterTraceByMixtureElementForCodonR)
		.method("getLambdaPrimeParameterTraceByMixtureElementForCodon", &RFPTrace::getLambdaPrimeParameterTraceByMixtureElementForCodonR)
		;

	class_<FONSETrace>("FONSETrace")
		.derives<Trace>("Trace")
		//These methods have specific R wrappers
		.method("getMutationParameterTraceByMixtureElementForCodon", &FONSETrace::getMutationParameterTraceByMixtureElementForCodonR)
		.method("getSelectionParameterTraceByMixtureElementForCodon", &FONSETrace::getSelectionParameterTraceByMixtureElementForCodonR)
		;
}
#endif

