
# Assumes values have been set by parent makefile:
#
# BOOST_ROOT : currently tested on Boost v1.44
#
# REDIST_DIR : location of other dependencies (besides boost)
#
# SAMTOOLS_ID : uses an older version of the samtools library (0.1.8) for bam parsing. This is included in redist/ and should not need to be set.
# TABIX_ID :
# CODEMIN_ID : 
#

IDS := REDIST_DIR BOOST_ID CODEMIN_ID SAMTOOLS_ID TABIX_ID

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

file-exists = $(wildcard $1)
check-exists-indirect = $(if $(call file-exists,$($1)),,$(error Variable '$1' refers to a non-existent path '$($1)'))
$(foreach I,REDIST_DIR, $(call check-exists-indirect,$I))

SAMTOOLS_DIR := $(REDIST_DIR)/$(SAMTOOLS_ID)
TABIX_DIR := $(REDIST_DIR)/$(TABIX_ID)
CODEMIN_INCLUDE_DIR := $(REDIST_DIR)/$(CODEMIN_ID)/include

WARNFLAGS := -Wall #-Werror
OPTFLAGS := -O2
DEBUGFLAGS := -ggdb -O0
XFLAGS := $(WARNFLAGS) $(OPTFLAGS)

BOOST_REDIST_ROOT := $(REDIST_DIR)/$(BOOST_ID)/stage
BOOST_INCLUDE_DIR := $(BOOST_REDIST_ROOT)/include
BOOST_LIB_DIR := $(BOOST_REDIST_ROOT)/lib

STARKA_INCLUDE_DIR:=$(CURDIR)/lib

CXXFLAGS := $(XFLAGS)
LDFLAGS := $(XFLAGS)

STAR_BIN := bin/starling2

VAR_BINS := $(STAR_BIN)

SHARED_CXX_DIR := blt_util blt_common common starling_common
SHARED_CXX_SRC := $(wildcard $(patsubst %,lib/%/*.cpp,$(SHARED_CXX_DIR)))
SHARED_OBJS := $(patsubst %.cpp,%.o,$(SHARED_CXX_SRC))
STAR_CXX_SRC := bin/starling2.cpp $(wildcard lib/starling/*.cpp) $(SHARED_CXX_SRC)
STAR_OBJS := $(patsubst %.cpp,%.o,$(STAR_CXX_SRC))

all: build

build: $(VAR_BINS)

install: build
	cp $(VAR_BINS) $(BIN_DIR)

$(VAR_BINS): CXXFLAGS += -I$(STARKA_INCLUDE_DIR) -I$(BOOST_INCLUDE_DIR) -I$(SAMTOOLS_DIR) -I$(TABIX_DIR) -I$(CODEMIN_INCLUDE_DIR)
$(VAR_BINS): LIBS := -lbam -ltabix -lz -lboost_program_options
$(VAR_BINS): LDLIBS += -L$(BOOST_LIB_DIR) -L$(SAMTOOLS_DIR) -L$(TABIX_DIR) $(LIBS)
$(VAR_BINS): CC := $(CXX)

$(STAR_BIN): $(STAR_OBJS) $(SHARED_OBJS)


clean:
	$(RM) $(STAR_OBJS) $(SHARED_OBJS) $(STAR_BIN)


# extra dev tags:
#
etags: TAGS

TAGS:
	find . -type f -iname "*.[ch]" -or -iname "*.cpp" -or -iname "*.hh" | sed "s/\.\///" | etags - 


