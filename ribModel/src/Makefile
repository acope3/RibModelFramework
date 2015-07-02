
roc := $(true)
fonse := $(false)
rfp := $(false)
 
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -fopenmp -D LINUX
endif
ifeq ($(UNAME_S),Darwin)
    CCFLAGS += -D OSX
endif
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
    CCFLAGS += -D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
    CCFLAGS += -D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
    CCFLAGS += -D ARM
endif

CXX ?= g++

SRC_PATH = .

SRC_EXT = cpp

COMPILE_FLAGS = -std=c++11 -Wall -Wextra -O2

INCLUDES = -I $(SRC_PATH)/include/

if roc
	INCLUDES += -I $(SRC_PATH)/include/ROC/
	
endif
