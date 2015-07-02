
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

BASE_SOURCES= src/base/Trace.cpp src/base/Model.cpp src/base/Parameter.cpp
ROC_SOURCES= src/ROC/ROCTrace.cpp src/ROC/ROCModel.cpp src/ROC/ROCParameter.cpp
#RFP_SOURCES= src/RFP/RFPTrace.cpp src/RFP/RFPModel.cpp src/RFP/RFPParameter.cpp
#FONSE_SOURCES= src/FONSE/FONSETrace.cpp src/FONSE/FONSEModel.cpp src/FONSE/FONSEParameter.cpp
RCPP_SOURCES= src/RCPP/RCPP_Model.cpp src/RCPP/RCPP_Parameter.cpp src/RCPP/RCPP_Trace.cpp
ROOT_SOURCES= src/CovarianceMatrix.cpp src/Gene.cpp src/Genome.cpp src/MCMCAlgorithm.cpp src/SequenceSummary.cpp 

SOURCES= $(BASE_SOURCES) $(ROC_SOURCES) $(RCPP_SOURCES) $(ROOT_SOURCES)

BASE_INCLUDES= -I./src/include/base/
ROC_INCLUDES= -I./src/include/ROC/
#RFP_INCLUDES= -I./include/RFP/
#FONSE_INCLUDES= -I./include/FONSE/
RCPP_INCLUDES= -I/usr/lib/R/site-library/Rcpp/include/ -I./include/RCPP/
ROOT_INCLUDES= -I./src/include/

INCLUDES= $(BASE_INCLUDES) $(ROC_INCLUDES) $(RCPP_INCLUDES) $(ROOT_INCLUDES)

ifeq ($(OS),Windows_NT)
	CC=cl
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CC=g++
		PKG_CPPFLAGS=-g -fopenmp -Wall -Wextra $(INCLUDES) -I../inst/include/ -DCPPTOML_USE_MAP
		PKG_LIBS= -fopenmp
	else
	ifeq ($(UNAME_S),Darwin)
		CC=clang++
		PKG_CPPFLAGS= -g -Wall -Wextra $(INCLUDES) -I../inst/include/ -DCPPTOML_USE_MAP
		PKG_LIBS=
	endif
	endif
endif

R_FILES := $(wildcard R/*.R)
SRC_FILES := $(SOURCES)
PKG_FILES := DESCRIPTION NAMESPACE $(R_FILES) $(SRC_FILES)

.PHONY: tarball install check clean build

tarball: $(PKG_NAME)_$(PKG_VERSION).tar.gz

$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	R CMD build .

check: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check $(PKG_NAME)_$(PKG_VERSION).tar.gz

build: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

clean: 
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -rf $(PKG_NAME).Rcheck
	-rm -rf man/*

.PHONY: list
list:
	@echo "R files:"
	@echo $(R_FILES)
	@echo "Source files:"
	@echo $(SRC_FILES)
