OS := $(shell uname)
ifeq ($(OS), Darwin)
    # mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib
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


.PHONY: default
default: fastcat bamstats bamindex


htslib/libhts.a:
	@echo Compiling $(@F)
	cd htslib/ \
		&& autoheader \
		&& autoconf \
		&& autoreconf --install \
		&& CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& make -j 4

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

fastcat: src/fastcat/main.o src/fastcat/args.o src/fastcat/writer.o src/fastqcomments.o src/common.o src/stats.o src/kh_counter.o $(STATIC_HTSLIB) zlib-ng/libz.a
	$(CC) -Isrc -Izlib-ng $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm $(EXTRA_LIBS) \
		-o $@

bamstats: src/bamstats/main.o src/bamstats/args.o src/bamstats/readstats.o src/bamstats/bamiter.o src/fastqcomments.o src/common.o src/stats.o src/kh_counter.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamindex: src/bamindex/main.o src/bamindex/build_main.o src/bamindex/fetch_main.o src/bamindex/dump_main.o src/bamindex/index.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

test/rg_parse: test/rg_parse.o src/common.o 
	$(CC) -Isrc $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm $(EXTRA_LIBS) \
		-o $@

.PHONY: clean
clean:
	rm -rf fastcat bamstats bamindex src/fastcat/*.o src/bamstats/*.o src/bamindex/*.o src/*.o

.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean

.PHONY: mem_check_fastcat
mem_check_fastcat: fastcat
	rm -rf fastcat-histograms
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./fastcat test/data/*.fastq.gz > /dev/null

.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: fastcat
	rm -rf demultiplex
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./fastcat test/data/*.fastq.gz --demultiplex demultiplex > /dev/null

.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: fastcat
	rm -rf demultiplex
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./fastcat test/data/*.fastq.gz --demultiplex demultiplex > /dev/null

.PHONY: mem_check_bamstats
mem_check_bamstats: bamstats
	rm -rf bamstats-histograms
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamstats test/bamstats/400ecoli-with-qcfail.bam
	rm -rf bamstats-histograms
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamstats test/parse_rg/dna_r10.4.1_e8.2_400bps_hac@v4.3.0.bam

.PHONY: mem_check_bamstats_duplex
mem_check_bamstats_duplex: bamstats
	[ -d bamstats-histograms ] && rm -rf bamstats-histograms
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamstats test/bamstats/310dx.bam

.PHONY: mem_check_bamindex-build
mem_check_bamindex-build: bamindex
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamindex build test/bamindex/400.bam

test/bamindex/400.bam.bci: bamindex
	./bamindex dump test/bamindex/400.bam.bci

.PHONY: mem_check_bamindex-dump
mem_check_bamindex-dump: bamindex test/bamindex/400.bam.bci
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamindex dump test/bamindex/400.bam.bci

.PHONY: mem_check_bamindex-fetch
mem_check_bamindex-fetch: bamindex test/bamindex/400.bam.bci
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -s \
		./bamindex fetch test/bamindex/400.bam --chunk 5 > /dev/null

.PHONY: mem_check_bamindex
mem_check_bamindex: mem_check_bamindex-build mem_check_bamindex-dump mem_check_bamindex-fetch

.PHONY: test_bamstats_NM
test_bamstats_NM: bamstats
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	../../bamstats ../bamstats_badNM/test.sam 2> err || grep "appears to contain implausible alignment information" err && rm -rf bamstats-histograms && \
	../../bamstats ../bamstats_zeroNM/test.sam
	rm -r test/test-tmp

.PHONY: test_bamstats_polya
test_bamstats_polya: bamstats
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	../../bamstats ../bamstats/RCS-100A.bam --poly_a > /dev/null && \
	diff bamstats-histograms/polya.hist ../bamstats/RCS-100A.bam.polya.hist
	rm -r test/test-tmp

.PHONY: regression_test_fastcat
regression_test_fastcat: fastcat
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	../../fastcat ../data -s sample --reheader -f per-file-stats.tsv -r per-read-stats.tsv \
		> concat.sorted.fastq && \
	bash -c 'diff <(sort per-file-stats.tsv) \
		<(sort ../fastcat_expected_results/per-file-stats.tsv)' && \
	bash -c 'diff <(sort per-read-stats.tsv) \
		<(sort ../fastcat_expected_results/per-read-stats.tsv)' && \
	bash -c "diff \
		<(cat concat.sorted.fastq | paste -d '|' - - - - | sort | tr '|' '\n') \
		<(${ZCAT} ../fastcat_expected_results/concat.reheader.sorted.fastq.gz | \
			paste -d '|' - - - - | sort | tr '|' '\n')" && \
	rm -rf fastcat-histograms/ && \
	../../fastcat ../data -s sample -f per-file-stats.tsv -r per-read-stats.tsv \
		> concat.sorted.fastq && \
	bash -c "diff \
		<(cat concat.sorted.fastq | paste -d '|' - - - - | sort | tr '|' '\n') \
		<(${ZCAT} ../fastcat_expected_results/concat.sorted.fastq.gz | \
			paste -d '|' - - - - | sort | tr '|' '\n')"
	rm -r test/test-tmp

# samtools fastq -T'*' sam2fastq/wf_basecalling_demo.sam > sam2fastq/wf_basecalling_demo.fastq
.PHONY: regression_test_sam_to_fastcat
regression_test_sam_to_fastcat: fastcat
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	../../fastcat ../sam2fastq/wf_basecalling_demo.fastq -s sample -H \
		--histograms fastcat-histograms1 -r per-read-stats.fastcat_once.tsv -f per-file-stats.fastcat_once.tsv \
		> wf_basecalling_demo.fastcat_once.fastq && \
	../../fastcat wf_basecalling_demo.fastcat_once.fastq -s sample -H \
		--histograms fastcat-histograms2 -r per-read-stats.fastcat_twice.tsv -f per-file-stats.fastcat_twice.tsv \
		> wf_basecalling_demo.fastcat_twice.fastq && \
	diff wf_basecalling_demo.fastcat_once.fastq wf_basecalling_demo.fastcat_twice.fastq && \
	bash -c 'diff <(sort -d per-file-stats.fastcat_once.tsv | sed 's,../sam2fastq/wf_basecalling_demo.fastq,FILE,') \
		<(sort -d per-file-stats.fastcat_twice.tsv | sed 's,wf_basecalling_demo.fastcat_once.fastq,FILE,')' && \
	bash -c 'diff <(sort per-read-stats.fastcat_once.tsv | sed 's,../sam2fastq/wf_basecalling_demo.fastq,FILE,') \
		<(sort per-read-stats.fastcat_twice.tsv | sed 's,wf_basecalling_demo.fastcat_once.fastq,FILE,')'
	rm -r test/test-tmp

.PHONY: regression_test_parse_rg_fastq
regression_test_parse_rg_fastq: fastcat
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	for i in ../parse_rg/*.fastq.gz; do \
		echo $$i; \
		../../fastcat $$i --histograms hist -l rg \
			> /dev/null; \
		diff rg $$i.callers; \
		rm -rf hist rg; \
	done;
	rm -r test/test-tmp

.PHONY: regression_test_parse_rg_bam
regression_test_parse_rg_bam: bamstats
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	for i in ../parse_rg/*.bam; do \
		../../bamstats $$i --histograms hist -l rg \
			> /dev/null; \
		diff rg $$i.callers; \
		rm -rf hist rg; \
	done;
	rm -r test/test-tmp

.PHONY: regression_test_parse_rd_fastq
regression_test_parse_rd_fastq: fastcat
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	for i in ../parse_rd/*.fastq; do \
		echo $$i; \
		../../fastcat $$i --histograms hist -i rd > /dev/null; \
		diff rd $$i.runids; \
		rm -rf hist rg; \
	done;
	rm -r test/test-tmp

.PHONY: regression_test_rg_parsing
regression_test_rg_parsing: test/rg_parse
	./test/rg_parse
