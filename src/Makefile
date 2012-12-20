
# Assumes values have been set by parent makefile:
#

CXX_DIR = $(CURDIR)/c++
PERL_DIR = $(CURDIR)/perl


all: build
	
build:
	$(MAKE) -C $(CXX_DIR) $@
	$(MAKE) -C $(PERL_DIR) $@

install:
	$(MAKE) -C $(CXX_DIR) $@
	$(MAKE) -C $(PERL_DIR) $@

clean:
	$(MAKE) -C $(CXX_DIR) $@
	$(MAKE) -C $(PERL_DIR) $@

