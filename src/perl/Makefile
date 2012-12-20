#
# Assumes values have been set by parent makefile:
#
IDS := INSTALL_DIR

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

file-exists = $(wildcard $1)
check-exists-indirect = $(if $(call file-exists,$($1)),,$(error Variable '$1' refers to a non-existent path '$($1)'))


all: build

check_install_dir:
	$(call check-exists-indirect,INSTALL_DIR)

install: check_install_dir
	rsync -a bin/ $(INSTALL_DIR)/bin
	rsync -a lib/ $(INSTALL_DIR)/lib
	rsync -a libexec/ $(INSTALL_DIR)/libexec

build:

clean:

