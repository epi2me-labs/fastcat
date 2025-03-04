OS := $(shell uname)
ifeq ($(OS), Darwin)
    # mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib -L$(shell brew --prefix curl)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
    ARGP_INC ?= -I$(shell brew --prefix argp-standalone)/include
    CFLAGS ?= -fpic -O3 ${ARGP_INC}
    ZCAT = "gzcat"
else
    ARGP ?=
    ARGP_INC ?=
    CFLAGS ?= -fpic -msse3 -O3 ${ARGP_INC}
    ZCAT = "zcat"
endif

VALGRIND ?= valgrind

CC ?= gcc
STATIC_HTSLIB ?= htslib/libhts.a
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
    CFLAGS += -DNOTHREADS
endif

# we can't do pedantic because min/max macros lead to:
#     "ISO C forbids braced-groups within expressions [-Werror=pedantic]"
ifeq ($(shell $(CC) --version | grep clang | wc -l), 0)
    WARNINGS = -Werror -Wall -Wextra -Wno-incompatible-pointer-types
else
    WARNINGS = -Werror -Wall -Wextra -Wpedantic -Wno-language-extension-token -Wno-gnu-statement-expression -Wno-incompatible-function-pointer-types
endif

GRIND = $(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s
ifeq ($(OS), Darwin)
	GRIND =
endif
# optionally run all tests under valgrind
ifeq ($(PEPPER), 1)
	PEPPER = $(GRIND)
else
	PEPPER = 
endif


.PHONY:
default: fastcat bamstats bamindex

.PHONY:
test: test_fastcat test_bamstats test_meta test_bamindex

.PHONY:
test_memory: mem_check_fastcat mem_check_bamstats mem_check_bamindex

.PHONY:
clean:
	rm -rf fastcat bamstats bamindex src/fastcat/*.o src/bamstats/*.o src/bamindex/*.o src/*.o

.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean


###
# build stages

htslib/libhts.a:
	@echo Compiling $(@F)
	cd htslib/ \
		&& autoheader \
		&& autoconf \
		&& autoreconf --install \
		&& CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& make -j 4

# just for testing
SAMVER=1.21
samtools:
	curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2;
	tar -xjf samtools-${SAMVER}.tar.bz2;
	rm samtools-${SAMVER}.tar.bz2
	cd samtools-${SAMVER} && make -j 4
	cp samtools-${SAMVER}/samtools $@

#TODO: for conda we could use zlib-ng from conda-forge

zlib-ng/zlib.h:
	@echo Configuring zlib-ng
	cd zlib-ng/ \
		&& CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure --zlib-compat \

zlib-ng/libz.a: zlib-ng/zlib.h
	@echo Compiling $(@F)
	cd zlib-ng/ \
		&& make -j 4 libz.a

src/%.o: src/%.c zlib-ng/zlib.h
	$(CC) -Isrc -Ihtslib -Izlib-ng -c -pthread $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $< -o $@

fastcat: src/version.o src/fastcat/main.o src/fastcat/args.o src/fastcat/writer.o src/fastqcomments.o src/common.o src/stats.o src/kh_counter.o $(STATIC_HTSLIB) zlib-ng/libz.a
	$(CC) -Isrc -Izlib-ng $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamstats: src/version.o src/bamstats/main.o src/bamstats/args.o src/bamstats/readstats.o src/bamstats/bamiter.o src/fastqcomments.o src/common.o src/regiter.o src/stats.o src/kh_counter.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamindex: src/version.o src/bamindex/main.o src/bamindex/build_main.o src/bamindex/fetch_main.o src/bamindex/dump_main.o src/bamindex/index.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

test/rg_parse: src/version.o test/rg_parse.o src/common.o 
	$(CC) -Isrc $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm $(EXTRA_LIBS) \
		-o $@


###
# fastcat tests

.PHONY:
test_fastcat: mem_check_fastcat mem_check_fastcat_demultiplex mem_check_fastcat_bam mem_check_fastcat_demultiplex_bam test_fastcat_bam_equivalent

.PHONY: mem_check_fastcat
mem_check_fastcat: fastcat
	rm -rf fastcat-histograms
	$(GRIND) ./fastcat test/data/*.fastq.gz > /dev/null

.PHONY: mem_check_fastcat_bam
mem_check_fastcat_bam: fastcat
	rm -rf fastcat-histograms
	$(GRIND) ./fastcat test/data/*.fastq.gz -B > /dev/null

.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: fastcat
	rm -rf demultiplex
	$(GRIND) ./fastcat test/data/*.fastq.gz --demultiplex demultiplex > /dev/null

.PHONY: mem_check_fastcat_demultiplex_bam
mem_check_fastcat_demultiplex_bam: fastcat
	rm -rf demultiplex
	$(GRIND) ./fastcat test/data/*.fastq.gz --demultiplex demultiplex -B > /dev/null

.PHONY: test_fastcat_bam_equivalent
fastcat_bam_equivalent: fastcat bamstats samtools
	@echo ""
	@echo "Testing fastcat bam equivalence"
	rm -rf test/test-tmp-fcb-equiv-van*
	rm -rf test/test-tmp-fcb-equiv-bam*
	$(PEPPER) ./fastcat test/data/*.fastq.gz --histograms test/test-tmp-fcb-equiv-van --reheader | ./samtools import -T '*' - | ./test/sort-sam.py > test/test-tmp-fcb-equiv-van.sam && \
	$(PEPPER) ./fastcat test/data/*.fastq.gz --histograms test/test-tmp-fcb-equiv-bam -B | ./samtools view | ./test/sort-sam.py > test/test-tmp-fcb-equiv-bam.sam && \
	diff test/test-tmp-fcb-equiv-van.sam test/test-tmp-fcb-equiv-bam.sam


###
# bamstats tests

.PHONY: 
test_bamstats: test_bamstats_NM test_bamstats_polya mem_check_bamstats

.PHONY: test_bamstats_NM
test_bamstats_NM: bamstats
	rm -rf test/test-tmp-bs-nm
	mkdir test/test-tmp-bs-nm && \
	cd test/test-tmp-bs-nm && \
	$(PEPPER) ../../bamstats ../bamstats_badNM/test.sam 2> err || grep "appears to contain implausible alignment information" err && rm -rf bamstats-histograms-bs-nm && \
	rm -rf bamstats-histograms && \
	$(PEPPER) ../../bamstats ../bamstats_zeroNM/test.sam
	rm -r test/test-tmp-bs-nm

.PHONY: test_bamstats_polya
test_bamstats_polya: bamstats
	rm -rf test/test-tmp-bs-pa
	mkdir test/test-tmp-bs-pa && \
	cd test/test-tmp-bs-pa && \
	$(PEPPER) ../../bamstats ../bamstats/RCS-100A.bam --poly_a > /dev/null && \
	diff bamstats-histograms/polya.hist ../bamstats/RCS-100A.bam.polya.hist
	rm -r test/test-tmp-bs-pa

.PHONY:
mem_check_bamstats: bamstats
	@echo "Memcheck bamstats with good data"
	rm -rf bamstats-histograms
	$(GRIND) ./bamstats test/parse_rg/dna_r10.4.1_e8.2_400bps_hac@v4.3.0.bam > /dev/null
	@echo "Memcheck bamstats with bad data"
	@echo ""
	rm -rf bamstats-histograms
	$(GRIND) ./bamstats test/parse_rg/bad-ones.bam > /dev/null
	@echo ""
	@echo "Memcheck bamstats with qcfails"
	rm -rf bamstats-histograms
	$(GRIND) ./bamstats test/bamstats/400ecoli-with-qcfail.bam > /dev/null
	@echo ""
	@echo "Memcheck bamstats duplex"
	rm -rf bamstats-histograms
	$(GRIND) ./bamstats test/bamstats/310dx.bam


###
# meta data tests (both fastcat and bamstats)

.PHONY:
test_meta: test_meta_fastcat test_meta_bamstats

.PHONY: test_meta_fastcat
test_meta_fastcat: fastcat
	rm -rf test/test-tmp-meta-fastq
	mkdir test/test-tmp-meta-fastq && \
	cd test/test-tmp-meta-fastq && \
	set -e; \
	for i in ../parse_rd/*.fastq; do \
		echo $$i; \
		$(PEPPER) ../../fastcat $$i --histograms hist -i rd > /dev/null; \
		diff rd $$i.runids || exit 1; \
		rm -rf hist rg; \
	done;
	rm -r test/test-tmp-meta-fastq

.PHONY: test_meta_bamstats
test_meta_bamstats: bamstats
	rm -rf test/test-tmp-meta-bam
	mkdir test/test-tmp-meta-bam && \
	cd test/test-tmp-meta-bam && \
	set -e; \
	for i in ../parse_rg/*.bam; do \
		$(PEPPER) ../../bamstats $$i --histograms hist -l rg \
			> /dev/null; \
		diff rg $$i.callers || exit 1; \
		rm -rf hist rg; \
	done;
	rm -r test/test-tmp-meta-bam

.PHONY: regression_test_rg_parsing
regression_test_rg_parsing: test/rg_parse
	$(PEPPER) ./test/rg_parse


###
# bamindex tests

.PHONY:
test_bamindex: mem_check_bamindex-build mem_check_bamindex-dump mem_check_bamindex-fetch 

.PHONY: mem_check_bamindex-build
mem_check_bamindex-build: bamindex
	$(GRIND) ./bamindex build test/bamindex/400.bam

.PHONY: mem_check_bamindex-dump
mem_check_bamindex-dump: bamindex mem_check_bamindex-build
	$(GRIND) ./bamindex dump test/bamindex/400.bam.bci > /dev/null

.PHONY: mem_check_bamindex-fetch
mem_check_bamindex-fetch: bamindex mem_check_bamindex-build
	$(GRIND) ./bamindex fetch test/bamindex/400.bam --chunk 5 > /dev/null

