#
# relies on 'IDS' definitions from parent makefile
#

IDS := BGZF_CAT_ID

BGZF_CAT_ID := bgzf_extras-v0.1

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

id-to-dir = $($1)
DIRS := $(foreach I,$(IDS), $(call id-to-dir,$I))

LINKS := bgzf_extras

.PHONY: $(IDS)

all: build
	
build: $(LINKS)

$(IDS): $(DIRS)

safe-expand = \
mkdir -p $1.extract.tmp; \
$2 -dc < $1 | tar -xf - -C $1.extract.tmp; \
mv $1.extract.tmp/* .; \
rmdir $1.extract.tmp

maybe-expand = $(if $(wildcard $1),$(call safe-expand,$1,$2,$3),)

$(DIRS):
	$(call maybe-expand,$@.tar.bz2,bzip2,$@)
	$(call maybe-expand,$@.tar.gz,gzip,$@)


# vcftools has a makefile issue preventing recursive make on some linux distros unless MAKEFLAGS is cleared first:
bgzf_extras: $(BGZF_CAT_ID)
	$(MAKE) -C $< && \
	ln -sf $< $@ 

clean:
	$(RM) -rf $(DIRS) $(LINKS)

