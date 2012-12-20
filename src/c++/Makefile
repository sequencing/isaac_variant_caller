#
# Assumes values have been set by parent makefile:
#
IDS := INSTALL_DIR

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

file-exists = $(wildcard $1)
check-exists-indirect = $(if $(call file-exists,$($1)),,$(error Variable '$1' refers to a non-existent path '$($1)'))


WARNFLAGS := -Wall
OPTFLAGS := -O2
XFLAGS := $(WARNFLAGS) $(OPTFLAGS)

CXXFLAGS := $(XFLAGS)
LDFLAGS := $(XFLAGS)

FASTA_BIN := libexec/countFastaBases

all: build

check_install_dir:
	$(call check-exists-indirect,INSTALL_DIR)


install: check_install_dir build
	cp $(FASTA_BIN) $(INSTALL_DIR)/$(FASTA_BIN)

build: $(FASTA_BIN)

$(FASTA_BIN): $(FASTA_BIN).cpp

clean:
	$(RM) $(FASTA_BIN) libexec/*.o

