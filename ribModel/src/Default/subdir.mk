################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CovarianceMatrix.cpp \
../Gene.cpp \
../Genome.cpp \
../MCMCAlgorithm.cpp \
../Model.cpp \
../Parameter.cpp \
../RCPP_Model.cpp \
../RCPP_Parameter.cpp \
../RCPP_Trace.cpp \
../RFPModel.cpp \
../RFPParameter.cpp \
../RFPTrace.cpp \
../ROCModel.cpp \
../ROCParameter.cpp \
../ROCTrace.cpp \
../SequenceSummary.cpp \
../Trace.cpp \
../main.cpp 

OBJS += \
./CovarianceMatrix.o \
./Gene.o \
./Genome.o \
./MCMCAlgorithm.o \
./Model.o \
./Parameter.o \
./RCPP_Model.o \
./RCPP_Parameter.o \
./RCPP_Trace.o \
./RFPModel.o \
./RFPParameter.o \
./RFPTrace.o \
./ROCModel.o \
./ROCParameter.o \
./ROCTrace.o \
./SequenceSummary.o \
./Trace.o \
./main.o 

CPP_DEPS += \
./CovarianceMatrix.d \
./Gene.d \
./Genome.d \
./MCMCAlgorithm.d \
./Model.d \
./Parameter.d \
./RCPP_Model.d \
./RCPP_Parameter.d \
./RCPP_Trace.d \
./RFPModel.d \
./RFPParameter.d \
./RFPTrace.d \
./ROCModel.d \
./ROCParameter.d \
./ROCTrace.d \
./SequenceSummary.d \
./Trace.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Users/roxasoath1/Desktop/RibModelFramework/ribModel/src/include -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


