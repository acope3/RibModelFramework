################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../FONSE/FONSEParameter.cpp \
../FONSE/FONSETrace.cpp 

OBJS += \
./FONSE/FONSEParameter.o \
./FONSE/FONSETrace.o 

CPP_DEPS += \
./FONSE/FONSEParameter.d \
./FONSE/FONSETrace.d 


# Each subdirectory must supply rules for building sources it contributes
FONSE/%.o: ../FONSE/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Users/roxasoath1/Desktop/RibModelFramework/ribModel/src/include -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


