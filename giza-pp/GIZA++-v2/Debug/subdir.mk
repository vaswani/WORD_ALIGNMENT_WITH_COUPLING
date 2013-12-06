################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ATables.cpp \
../AlignTables.cpp \
../Dictionary.cpp \
../ForwardBackward.cpp \
../HMMTables.cpp \
../MoveSwapMatrix.cpp \
../NTables.cpp \
../Parameter.cpp \
../Perplexity.cpp \
../TTables.cpp \
../alignment.cpp \
../collCounts.cpp \
../getSentence.cpp \
../hmm.cpp \
../logprob.cpp \
../main.cpp \
../model1.cpp \
../model2.cpp \
../model2to3.cpp \
../model3.cpp \
../model345-peg.cpp \
../model3_viterbi.cpp \
../model3_viterbi_with_tricks.cpp \
../myassert.cpp \
../parse.cpp \
../plain2snt.cpp \
../reports.cpp \
../small_snt2cooc.cpp \
../snt2cooc.cpp \
../snt2plain.cpp \
../transpair_model3.cpp \
../transpair_model4.cpp \
../transpair_model5.cpp \
../utility.cpp \
../vocab.cpp 

OBJS += \
./ATables.o \
./AlignTables.o \
./Dictionary.o \
./ForwardBackward.o \
./HMMTables.o \
./MoveSwapMatrix.o \
./NTables.o \
./Parameter.o \
./Perplexity.o \
./TTables.o \
./alignment.o \
./collCounts.o \
./getSentence.o \
./hmm.o \
./logprob.o \
./main.o \
./model1.o \
./model2.o \
./model2to3.o \
./model3.o \
./model345-peg.o \
./model3_viterbi.o \
./model3_viterbi_with_tricks.o \
./myassert.o \
./parse.o \
./plain2snt.o \
./reports.o \
./small_snt2cooc.o \
./snt2cooc.o \
./snt2plain.o \
./transpair_model3.o \
./transpair_model4.o \
./transpair_model5.o \
./utility.o \
./vocab.o 

CPP_DEPS += \
./ATables.d \
./AlignTables.d \
./Dictionary.d \
./ForwardBackward.d \
./HMMTables.d \
./MoveSwapMatrix.d \
./NTables.d \
./Parameter.d \
./Perplexity.d \
./TTables.d \
./alignment.d \
./collCounts.d \
./getSentence.d \
./hmm.d \
./logprob.d \
./main.d \
./model1.d \
./model2.d \
./model2to3.d \
./model3.d \
./model345-peg.d \
./model3_viterbi.d \
./model3_viterbi_with_tricks.d \
./myassert.d \
./parse.d \
./plain2snt.d \
./reports.d \
./small_snt2cooc.d \
./snt2cooc.d \
./snt2plain.d \
./transpair_model3.d \
./transpair_model4.d \
./transpair_model5.d \
./utility.d \
./vocab.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


