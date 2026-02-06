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
CFLAGS += -MMD -MP

VALGRIND ?= valgrind
DEFLATEVER=1.24

CC ?= gcc
AR ?= ar
STATIC_HTSLIB ?= htslib/libhts.a
HTSLIB_PATCH ?= patches/htslib-aux-parse.patch
HTSLIB_PATCH_STAMP ?= .fastcat_patch_applied
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
    CFLAGS += -DNOTHREADS
endif

LIBDEF=$(abspath libdeflate-$(DEFLATEVER))
LIBDEF_LIB=$(abspath libdeflate-$(DEFLATEVER)/build)
HTS_BUILD_FLAGS =
STATIC_LIBDEFLATE =

ifeq ($(USE_DEFLATE), 1)
    STATIC_HTSLIB += libdeflate-$(DEFLATEVER)/build/libdeflate.a
    STATIC_LIBDEFLATE = libdeflate-$(DEFLATEVER)/build/libdeflate.a
    HTS_BUILD_FLAGS = CPPFLAGS="-I$(LIBDEF)" LDFLAGS="-L$(LIBDEF_LIB)" LIBS="-ldeflate"
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

SHELL := /bin/bash
REPO_ROOT := $(CURDIR)
TEST_ROOT := $(REPO_ROOT)/test
TEST_EXPECTED_FASTCAT := $(TEST_ROOT)/fastcat_expected_results
TEST_EXPECTED_DEMULTIPLEX := $(TEST_EXPECTED_FASTCAT)/demultiplex
TEST_TMP_DIR ?= $(TEST_ROOT)/tmp
TEST_TMP_FASTCAT = $(TEST_TMP_DIR)/fastcat
TEST_TMP_BAMSTATS = $(TEST_TMP_DIR)/bamstats
TEST_TMP_META = $(TEST_TMP_DIR)/meta
TEST_TMP_BAMINDEX = $(TEST_TMP_DIR)/bamindex
TEST_TMP_FASTLINT = $(TEST_TMP_DIR)/fastlint
TEST_TMP_BAMCOVERAGE = $(TEST_TMP_DIR)/bamcoverage

comma := ,

RUN_IN_TMP = rm -rf $(1); mkdir -p $(1) && pushd $(1) > /dev/null && $(2) && popd > /dev/null && rm -rf $(1)
COMPARE_TABLES = \
	head -n 1 $(1) > $(3).header && \
	head -n 1 $(2) > $(3).expected.header && \
	diff $(3).header $(3).expected.header && \
	tail -n +2 $(1) | LC_ALL=C sort > $(3).sorted && \
	tail -n +2 $(2) | LC_ALL=C sort > $(3).expected.sorted && \
	diff $(3).sorted $(3).expected.sorted
COMPARE_TABLES_SKIP_FIRST = \
	head -n 1 $(1) | cut -f2- > $(3).header && \
	head -n 1 $(2) | cut -f2- > $(3).expected.header && \
	diff $(3).header $(3).expected.header && \
	tail -n +2 $(1) | cut -f2- | LC_ALL=C sort > $(3).sorted && \
	tail -n +2 $(2) | cut -f2- | LC_ALL=C sort > $(3).expected.sorted && \
	diff $(3).sorted $(3).expected.sorted
COMPARE_FASTQ = \
	case "$(1)" in *.gz) gunzip -c "$(1)" ;; *) cat "$(1)" ;; esac \
		| $(TEST_ROOT)/sort-fastq.py > $(3).sorted && \
	case "$(2)" in *.gz) gunzip -c "$(2)" ;; *) cat "$(2)" ;; esac \
		| $(TEST_ROOT)/sort-fastq.py > $(3).expected.sorted && \
	diff $(3).sorted $(3).expected.sorted
TEST_DESC = printf '\n%s\n' "$(1)"
RUN_FASTCAT = $(PEPPER) $(REPO_ROOT)/fastcat
RUN_BAMSTATS = $(PEPPER) $(REPO_ROOT)/bamstats


.PHONY:
default: fastcat bamstats bamindex fastlint

.PHONY:
test: test_unit test_fastcat test_bamstats test_meta test_bamindex test_fastlint test_bamcoverage

.PHONY:
test_memory: mem_check_fastcat mem_check_bamstats mem_check_bamindex mem_check_fastlint mem_check_bamcoverage

.PHONY:
clean:
	rm -rf fastcat bamstats bamindex libfastcatcore.a src/fastcat/*.o src/fastlint/*.o src/bamstats/*.o src/bamindex/*.o src/*.o src/sink/*.o

.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean

# core library
CORE_OBJS = src/common.o src/fc_log.o src/stats.o src/fc_stats_engine.o src/fastqcomments.o src/fc_read.o src/fc_filter.o src/fc_bam_utils.o src/fc_dust.o src/fc_iter.o src/fc_writer.o src/fastcat/parsing.o src/sink/common.o src/sink/output.o src/sink/perread.o src/sink/perfile.o src/sink/counter.o src/sink/hist.o src/sink/coverage.o src/sink/flagstats.o src/bamcoverage/coverage.o src/sdust/sdust.o src/sdust/kalloc.o src/kh_counter.o src/regiter.o
CORE_LIB = libfastcatcore.a

# unit tests
UNIT_TEST_SRCS = test/unit/main.c test/unit/test_common.c test/unit/test_stats.c test/unit/test_stats_engine.c test/unit/test_fastqcomments.c test/unit/test_filters.c test/unit/test_rg_parse.c test/unit/test_regiter.c test/unit/test_quality.c test/unit/test_bam_utils.c test/unit/test_qualcalcs.c test/unit/test_iter.c test/unity/unity.c
UNITY_DEFS = -DUNITY_INCLUDE_DOUBLE -DUNITY_DOUBLE_PRECISION=1e-6
UNIT_TEST_OBJS = $(UNIT_TEST_SRCS:.c=.o)

test/unit/%.o: test/unit/%.c
	$(CC) -Isrc -Ihtslib -Izlib-ng -Itest/unity $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@

test/unity/%.o: test/unity/%.c
	$(CC) -Itest/unity $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@

test/unit/unit_tests: $(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib -Izlib-ng $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

.PHONY: test_unit
test_unit: test/unit/unit_tests
	@$(call TEST_DESC,Unit tests)
	$(PEPPER) ./test/unit/unit_tests


###
# build stages



$(HTSLIB_PATCH_STAMP): $(HTSLIB_PATCH)
	@if [ -f "$(HTSLIB_PATCH)" ] && { ! grep -q "sam_aux_parse(" htslib/htslib/sam.h || ! grep -q "sam_format_aux(" htslib/htslib/sam.h; }; then \
		echo "Applying htslib patch $(HTSLIB_PATCH)"; \
		cd htslib && patch -p1 -N < ../$(HTSLIB_PATCH); \
	fi
	@touch $(HTSLIB_PATCH_STAMP)

htslib/libhts.a: $(STATIC_LIBDEFLATE) | $(HTSLIB_PATCH_STAMP)
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd htslib/ \
		&& autoheader \
		&& autoconf \
		&& autoreconf --install \
		&& $(HTS_BUILD_FLAGS)  CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& make -j 4


libdeflate-$(DEFLATEVER)/build/libdeflate.a:
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	curl -L -o libdeflate-v${DEFLATEVER}.tar.gz https://github.com/ebiggers/libdeflate/archive/refs/tags/v${DEFLATEVER}.tar.gz \
		&& tar -xzf libdeflate-v${DEFLATEVER}.tar.gz \
		&& cd libdeflate-${DEFLATEVER} \
		&& cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
			-DCMAKE_INSTALL_PREFIX="$(pwd)/install" \
		&& cmake --build build -j \


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

src/fastqcomments.o: | $(HTSLIB_PATCH_STAMP)
src/fc_read.o: | $(HTSLIB_PATCH_STAMP)

$(CORE_LIB): $(CORE_OBJS) Makefile
	rm -f $@
	$(AR) rcs $@ $(CORE_OBJS)

-include $(wildcard src/*.d)
-include $(wildcard src/sink/*.d)
-include $(wildcard test/unit/*.d)
-include $(wildcard test/unity/*.d)

fastcat: src/version.o src/fastcat/main.o src/fastcat/args.o $(CORE_LIB) $(STATIC_HTSLIB) zlib-ng/libz.a
	$(CC) -Isrc -Izlib-ng $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

fastlint: src/version.o src/fastlint/main.o src/fastlint/args.o src/sdust/sdust.o src/sdust/kalloc.o src/fc_log.o $(STATIC_HTSLIB) zlib-ng/libz.a
	$(CC) -Isrc -Izlib-ng $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamstats: src/version.o src/bamstats/main.o src/bamstats/args.o $(CORE_LIB) $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamcoverage: src/version.o src/bamcoverage/main.o src/bamcoverage/args.o src/bamcoverage/coverage.o src/common.o src/regiter.o src/fc_log.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamindex: src/version.o src/bamindex/main.o src/bamindex/build_main.o src/bamindex/fetch_main.o src/bamindex/dump_main.o src/bamindex/index.o src/fc_log.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib $(WARNINGS) -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

###
# fastcat tests

.PHONY:
test_fastcat: mem_check_fastcat mem_check_fastcat_demultiplex mem_check_fastcat_bam mem_check_fastcat_demultiplex_bam test_fastcat_bam_equivalent test_fastcat_reheader_roundtrip test_fastcat_filters test_fastcat_demultiplex_outputs test_fastcat_demultiplex_split test_fastcat_demux_single_equiv test_fastcat_regression_results test_fastcat_perfile_all_filtered

.PHONY: test_fastcat_demux_single_equiv
test_fastcat_demux_single_equiv: fastcat
	@$(call TEST_DESC,Fastcat demux vs non-demux (single barcode))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/demux_single, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output plain_out --per-read > /dev/null && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output demux_out --per-read --demultiplex > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/runids.tsv,demux_out/runids.tsv,demux_single_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/basecallers.tsv,demux_out/basecallers.tsv,demux_single_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/per_file.tsv,demux_out/per_file.tsv,demux_single_perfile) && \
		test $$(tail -n +2 plain_out/per_read.tsv | wc -l) -eq $$(tail -n +2 demux_out/barcode0001/per_read.tsv | wc -l))

.PHONY: mem_check_fastcat
mem_check_fastcat: fastcat
	@$(call TEST_DESC,Memcheck fastcat (fastq))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/fastq, \
		$(GRIND) $(REPO_ROOT)/fastcat $(TEST_ROOT)/data/*.fastq.gz > /dev/null)

.PHONY: mem_check_fastcat_bam
mem_check_fastcat_bam: fastcat
	@$(call TEST_DESC,Memcheck fastcat (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/bam, \
		$(GRIND) $(REPO_ROOT)/fastcat $(TEST_ROOT)/data/*.fastq.gz -B > /dev/null)

.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: fastcat
	@$(call TEST_DESC,Memcheck fastcat demultiplex (fastq))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/demux, \
		$(GRIND) $(REPO_ROOT)/fastcat $(TEST_ROOT)/data/*.fastq.gz --demultiplex --output demultiplex > /dev/null)

.PHONY: mem_check_fastcat_demultiplex_bam
mem_check_fastcat_demultiplex_bam: fastcat
	@$(call TEST_DESC,Memcheck fastcat demultiplex (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/demux_bam, \
		$(GRIND) $(REPO_ROOT)/fastcat $(TEST_ROOT)/data/*.fastq.gz --demultiplex --output demultiplex -B > /dev/null)

.PHONY: test_fastcat_bam_equivalent
test_fastcat_bam_equivalent: fastcat bamstats samtools
	@$(call TEST_DESC,Fastcat BAM equivalence)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/*.fastq.gz --output out_reheader --reheader -B | $(REPO_ROOT)/samtools view | $(TEST_ROOT)/sort-sam.py > reheader.sam && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/*.fastq.gz --output out_bam -B | $(REPO_ROOT)/samtools view | $(TEST_ROOT)/sort-sam.py > bam.sam && \
		diff reheader.sam bam.sam)

.PHONY: test_fastcat_reheader_roundtrip
test_fastcat_reheader_roundtrip: fastcat samtools
	@$(call TEST_DESC,Fastcat reheader roundtrip)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/reheader, \
		$(RUN_FASTCAT) $(TEST_ROOT)/reheader/kv.fastq --output out_hist --reheader > out.fastq && \
		$(RUN_FASTCAT) $(TEST_ROOT)/reheader/kv.fastq --output out_orig -B | $(REPO_ROOT)/samtools view > orig.sam && \
		$(RUN_FASTCAT) out.fastq --output out_reheader -B | $(REPO_ROOT)/samtools view > reheader.sam && \
		diff orig.sam reheader.sam)

.PHONY: test_fastcat_filters
test_fastcat_filters: fastcat
	@$(call TEST_DESC,Fastcat filters)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/filters, \
		$(RUN_FASTCAT) $(TEST_ROOT)/filtering/fastcat_filters.fastq -a 5 -b 8 -q 20 --output out_filters > out.fastq && \
		grep -q "^@read_ok" out.fastq && \
		! grep -q "^@read_short" out.fastq && \
		! grep -q "^@read_long" out.fastq && \
		! grep -q "^@read_lowq" out.fastq)

.PHONY: test_fastcat_demultiplex_outputs
test_fastcat_demultiplex_outputs: fastcat
	@$(call TEST_DESC,Fastcat demultiplex outputs)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/demultiplex/outputs, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc0.fastq.gz $(TEST_ROOT)/data/bc1.fastq.gz $(TEST_ROOT)/data/bc2.fastq.gz \
			$(TEST_ROOT)/data/bcEmpty.fastq.gz $(TEST_ROOT)/data/bcMixed.fastq.gz --demultiplex --output out_demux > /dev/null && \
		$(call COMPARE_FASTQ,out_demux/unclassified/unclassified.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/unclassified.fastq.gz,demux_unclassified) && \
		$(call COMPARE_FASTQ,out_demux/barcode0001/barcode0001.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001.fastq.gz,demux_barcode0001) && \
		$(call COMPARE_FASTQ,out_demux/barcode0002/barcode0002.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0002.fastq.gz,demux_barcode0002) && \
		test -f out_demux/unclassified/unclassified.length.hist && \
		test -f out_demux/barcode0001/barcode0001.length.hist && \
		test -f out_demux/barcode0001/barcode0001.quality.hist)

.PHONY: test_fastcat_demultiplex_split
test_fastcat_demultiplex_split: fastcat
	@$(call TEST_DESC,Fastcat demultiplex split files)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/demultiplex/split, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --demultiplex --output out_split --reads-per-file 4 > /dev/null && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0000.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0000.fastq.gz,demux_split_0000) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0001.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0001.fastq.gz,demux_split_0001) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0002.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0002.fastq.gz,demux_split_0002) && \
		test -f out_split/barcode0001/barcode0001.length.hist)

.PHONY: test_fastcat_regression_results
test_fastcat_regression_results: fastcat
	@$(call TEST_DESC,Fastcat regression results)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/regression, \
		ln -s $(TEST_ROOT)/data data && mkdir -p work && \
		pushd work > /dev/null && \
		$(RUN_FASTCAT) ../data/*.fastq.gz --sample sample \
			--output $(TEST_TMP_FASTCAT)/regression/out1 \
			--per-read \
			> $(TEST_TMP_FASTCAT)/regression/concat.fastq && \
		$(RUN_FASTCAT) ../data/*.fastq.gz --reheader \
			--output $(TEST_TMP_FASTCAT)/regression/out2 \
			> $(TEST_TMP_FASTCAT)/regression/concat.reheader.fastq && \
		popd > /dev/null && \
		$(call COMPARE_TABLES,$(TEST_TMP_FASTCAT)/regression/out1/per_read.tsv,$(TEST_EXPECTED_FASTCAT)/per-read-stats.tsv,per-read) && \
		$(call COMPARE_TABLES,$(TEST_TMP_FASTCAT)/regression/out1/per_file.tsv,$(TEST_EXPECTED_FASTCAT)/per-file-stats.tsv,per-file) && \
		$(call COMPARE_FASTQ,concat.fastq,$(TEST_EXPECTED_FASTCAT)/concat.sorted.fastq.gz,concat) && \
		$(call COMPARE_FASTQ,concat.reheader.fastq,$(TEST_EXPECTED_FASTCAT)/concat.reheader.sorted.fastq.gz,concat.reheader))

.PHONY: test_fastcat_perfile_all_filtered
test_fastcat_perfile_all_filtered: fastcat
	@$(call TEST_DESC,Fastcat per-file stats when all reads filtered)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/perfile_filtered, \
		printf "@r1\nACGT\n+\n!!!!\n@r2\nACGT\n+\n!!!!\n" > short.fastq && \
		$(RUN_FASTCAT) short.fastq --min-length 100 --sample sample --output out > /dev/null && \
		awk 'NR==1{for(i=1;i<=NF;i++)h[$$i]=i;next} NR==2{ \
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_ok"]!=2||$$h["r_too_short"]!=2){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)

###
# bamstats tests

.PHONY: 
test_bamstats: test_bamstats_NM test_bamstats_polya test_bamstats_filters test_bamstats_vs_fastcat_unaligned test_bamstats_vs_fastcat_unaligned_demux test_bamstats_demux_single_equiv test_bamstats_perfile_all_filtered mem_check_bamstats
.PHONY: test_bamstats_coverage_demux_layout
test_bamstats: test_bamstats_coverage_demux_layout

.PHONY: test_bamstats_demux_single_equiv
test_bamstats_demux_single_equiv: fastcat bamstats
	@$(call TEST_DESC,Bamstats demux vs non-demux (single barcode))
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/demux_single, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output tmp_fc --per-read -B > bc1.bam && \
		$(RUN_BAMSTATS) bc1.bam --per-read --output plain_out > /dev/null && \
		$(RUN_BAMSTATS) bc1.bam --per-read --demultiplex --output demux_out > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/runids.tsv,demux_out/runids.tsv,bs_demux_single_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/basecallers.tsv,demux_out/basecallers.tsv,bs_demux_single_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,plain_out/per_file.tsv,demux_out/per_file.tsv,bs_demux_single_perfile) && \
		$(call COMPARE_TABLES,plain_out/per_read.tsv,demux_out/barcode0001/per_read.tsv,bs_demux_single_perread) && \
		diff plain_out/histograms/length.hist demux_out/barcode0001/barcode0001.length.hist && \
		diff plain_out/histograms/quality.hist demux_out/barcode0001/barcode0001.quality.hist)

.PHONY: test_bamstats_coverage_demux_layout
test_bamstats_coverage_demux_layout: bamstats samtools
	@$(call TEST_DESC,Bamstats coverage demux layout)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/coverage_demux, \
		cp $(TEST_ROOT)/bamstats/coverage_demux.sam cov.sam && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --coverage --demultiplex --output out > /dev/null && \
		if [ -d out/coverage ]; then echo "Unexpected top-level coverage directory"; exit 1; fi && \
		test -f out/barcode0001/coverage/global.bed.gz && \
		test -f out/barcode0002/coverage/global.bed.gz && \
		test ! -d out/unclassified/coverage)

.PHONY: test_bamstats_NM
test_bamstats_NM: bamstats
	@$(call TEST_DESC,Bamstats NM validation)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/nm, \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_badNM/test.sam --output out_bad 2> err || grep "appears to contain implausible alignment information" err && \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_zeroNM/test.sam --output out_zero)

.PHONY: test_bamstats_polya
test_bamstats_polya: bamstats
	@$(call TEST_DESC,Bamstats polya)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/polya, \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats/RCS-100A.bam --poly-a --output out_polya > /dev/null && \
		diff out_polya/histograms/polya.hist $(TEST_ROOT)/bamstats/RCS-100A.bam.polya.hist)

.PHONY: test_bamstats_filters
test_bamstats_filters: bamstats
	@$(call TEST_DESC,Bamstats filters)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/filters, \
		$(RUN_BAMSTATS) $(TEST_ROOT)/filtering/bamstats_filters.sam --tag-name HP --tag-value 1 --per-read --output out_tag > /dev/null && \
		grep -q "^read1" out_tag/per_read.tsv && \
		! grep -q "^read2" out_tag/per_read.tsv && \
		$(RUN_BAMSTATS) $(TEST_ROOT)/filtering/bamstats_filters.sam --read-group rg2 --per-read --output out_rg > /dev/null && \
		grep -q "^read2" out_rg/per_read.tsv && \
		! grep -q "^read1" out_rg/per_read.tsv && \
		$(RUN_BAMSTATS) $(TEST_ROOT)/filtering/bamstats_filters.sam --tag-name HP --tag-value 1 --read-group rg1 --per-read --output out_both > /dev/null && \
		grep -q "^read1" out_both/per_read.tsv && \
		! grep -q "^read2" out_both/per_read.tsv)

.PHONY: test_bamstats_vs_fastcat_unaligned
test_bamstats_vs_fastcat_unaligned: fastcat bamstats
	@$(call TEST_DESC,Bamstats vs fastcat (unaligned))
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/vs_fastcat, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bcMixed.fastq.gz --output fc_out --per-read -B > reads.bam && \
		$(RUN_BAMSTATS) reads.bam --per-read --output bs_out > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/runids.tsv,bs_out/runids.tsv,runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/basecallers.tsv,bs_out/basecallers.tsv,basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/per_file.tsv,bs_out/per_file.tsv,per_file) && \
		cut -f1$(comma)3-8 fc_out/per_read.tsv > pr.fc && \
		cut -f1$(comma)3-8 bs_out/per_read.tsv > pr.bs && \
		$(call COMPARE_TABLES,pr.fc,pr.bs,per_read_non_demux) && \
		diff fc_out/histograms/length.hist bs_out/histograms/length.unmap.hist && \
		diff fc_out/histograms/quality.hist bs_out/histograms/quality.unmap.hist)

.PHONY: test_bamstats_vs_fastcat_unaligned_demux
test_bamstats_vs_fastcat_unaligned_demux: fastcat bamstats
	@$(call TEST_DESC,Bamstats vs fastcat (unaligned, demux))
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/vs_fastcat_demux, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bcMixed.fastq.gz --output fc_out --per-read --demultiplex > /dev/null && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bcMixed.fastq.gz --per-read -B > reads.bam && \
		$(RUN_BAMSTATS) reads.bam --demultiplex --per-read --output bs_out > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/runids.tsv,bs_out/runids.tsv,runids_demux) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/basecallers.tsv,bs_out/basecallers.tsv,basecallers_demux) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/per_file.tsv,bs_out/per_file.tsv,per_file_demux) && \
		for f in $$(find fc_out -name per_read.tsv); do \
			rel=$${f#fc_out/}; \
			g=bs_out/$$rel; \
			if [ ! -f "$$g" ]; then echo "Missing $$g"; exit 1; fi; \
			cut -f1$(comma)3-8 "$$f" > pr.fc; \
			cut -f1$(comma)3-8 "$$g" > pr.bs; \
			$(call COMPARE_TABLES,pr.fc,pr.bs,per_read_$${rel//\//_}); \
			rm -f pr.fc pr.bs pr.header pr.expected.header pr.sorted pr.expected.sorted; \
		done && \
		for g in $$(find bs_out -name per_read.tsv); do \
			rel=$${g#bs_out/}; \
			f=fc_out/$$rel; \
			if [ ! -f "$$f" ]; then echo "Extra per_read.tsv in bamstats: $$rel"; exit 1; fi; \
		done && \
		for h in $$(find fc_out -name "*length.hist" -o -name "*quality.hist"); do \
			rel=$${h#fc_out/}; \
			base=$${rel%.hist}; \
			g=bs_out/$$base.unmap.hist; \
			if [ ! -f "$$g" ]; then echo "Missing $$g"; exit 1; fi; \
			diff "$$h" "$$g"; \
		done)

.PHONY: test_bamstats_perfile_all_filtered
test_bamstats_perfile_all_filtered: bamstats
	@$(call TEST_DESC,Bamstats per-file stats when all reads filtered)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/perfile_filtered, \
		printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n" > short.sam && \
		$(RUN_BAMSTATS) short.sam --tag-name ZZ --tag-value 1 --output out > /dev/null && \
		awk 'NR==1{for(i=1;i<=NF;i++)h[$$i]=i;next} NR==2{ \
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_ok"]!=2||$$h["r_too_short"]!=0){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)

.PHONY:
mem_check_bamstats: mem_check_bamstats_good mem_check_bamstats_bad mem_check_bamstats_qcfail mem_check_bamstats_duplex mem_check_bamstats_coverage

.PHONY: mem_check_bamstats_good
mem_check_bamstats_good: bamstats
	@$(call TEST_DESC,Memcheck bamstats with good data)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/good, \
		$(GRIND) $(REPO_ROOT)/bamstats $(TEST_ROOT)/parse_rg/dna_r10.4.1_e8.2_400bps_hac@v4.3.0.bam > /dev/null)

.PHONY: mem_check_bamstats_bad
mem_check_bamstats_bad: bamstats
	@$(call TEST_DESC,Memcheck bamstats with bad data)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/bad, \
		$(GRIND) $(REPO_ROOT)/bamstats $(TEST_ROOT)/parse_rg/bad-ones.bam > /dev/null)

.PHONY: mem_check_bamstats_qcfail
mem_check_bamstats_qcfail: bamstats
	@$(call TEST_DESC,Memcheck bamstats qcfail)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/qcfail, \
		$(GRIND) $(REPO_ROOT)/bamstats $(TEST_ROOT)/bamstats/400ecoli-with-qcfail.bam > /dev/null)

.PHONY: mem_check_bamstats_duplex
mem_check_bamstats_duplex: bamstats
	@$(call TEST_DESC,Memcheck bamstats duplex)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/duplex, \
		$(GRIND) $(REPO_ROOT)/bamstats $(TEST_ROOT)/bamstats/310dx.bam > /dev/null)

.PHONY: mem_check_bamstats_coverage
mem_check_bamstats_coverage: bamstats
	@$(call TEST_DESC,Memcheck bamstats coverage)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/coverage, \
		$(GRIND) $(REPO_ROOT)/bamstats --coverage $(TEST_ROOT)/bamstats/400ecoli.bam > /dev/null)

.PHONY:
bamstats_overwrite: bamstats
	@$(call TEST_DESC,Bamstats basecaller file exists)
	@$(PEPPER) cat test/bamstats/400ecoli.bam \
		| ./bamstats -u --output bamstats_results - \
		> /dev/null; \
		status=$$?; \
		if [ $$status -eq 1 ]; then \
			echo "Expected failure: bamstats basecaller file already exists"; \
		else \
			echo "Unexpected success: bamstats basecaller file already exists"; \
			exit 1; \
		fi

###
# meta data tests (both fastcat and bamstats)

.PHONY:
test_meta: test_meta_fastcat test_meta_bamstats

.PHONY: test_meta_fastcat
test_meta_fastcat: fastcat
	@$(call TEST_DESC,Metadata fastcat)
	$(call RUN_IN_TMP,$(TEST_TMP_META)/fastq, \
		ln -s $(TEST_ROOT)/parse_rd parse_rd && mkdir -p work && \
		pushd work > /dev/null && \
		set -e; \
		for i in ../parse_rd/*.fastq; do \
			echo $$i; \
			$(RUN_FASTCAT) $$i --output out --per-read > /dev/null; \
			diff out/runids.tsv $$i.runids || exit 1; \
			rm -rf out; \
		done; \
		popd > /dev/null)

.PHONY: test_meta_bamstats
test_meta_bamstats: bamstats
	@$(call TEST_DESC,Metadata bamstats)
	$(call RUN_IN_TMP,$(TEST_TMP_META)/bam, \
		ln -s $(TEST_ROOT)/parse_rg parse_rg && mkdir -p work && \
		pushd work > /dev/null && \
		set -e; \
		for i in ../parse_rg/*.bam; do \
			$(RUN_BAMSTATS) $$i --output out > /dev/null; \
			diff out/basecallers.tsv $$i.callers || exit 1; \
			rm -rf out; \
		done; \
		popd > /dev/null)

.PHONY: regression_test_rg_parsing
regression_test_rg_parsing: test/rg_parse
	@$(call TEST_DESC,RG parsing regression)
	$(PEPPER) ./test/rg_parse

###
# bamindex tests

.PHONY:
test_bamindex: mem_check_bamindex-build mem_check_bamindex-dump mem_check_bamindex-fetch 

.PHONY: mem_check_bamindex-build
mem_check_bamindex-build: bamindex
	@$(call TEST_DESC,Memcheck bamindex build)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/build, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/bamindex build 400.bam > /dev/null)

.PHONY: mem_check_bamindex-dump
mem_check_bamindex-dump: bamindex
	@$(call TEST_DESC,Memcheck bamindex dump)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/dump, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/bamindex build 400.bam > /dev/null && \
		$(GRIND) $(REPO_ROOT)/bamindex dump 400.bam.bci > /dev/null)

.PHONY: mem_check_bamindex-fetch
mem_check_bamindex-fetch: bamindex
	@$(call TEST_DESC,Memcheck bamindex fetch)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/fetch, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/bamindex build 400.bam > /dev/null && \
		$(GRIND) $(REPO_ROOT)/bamindex fetch 400.bam --chunk 5 > /dev/null)

###
# fastlint tests

test_fastlint: mem_check_fastlint

.PHONY:
mem_check_fastlint: fastlint
	@$(call TEST_DESC,Memcheck fastlint)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTLINT)/memcheck, \
		$(GRIND) $(REPO_ROOT)/fastlint $(TEST_ROOT)/data/*.fastq.gz > /dev/null)

###
# bamcoverage tests

test_bamcoverage: mem_check_bamcoverage

.PHONY:
mem_check_bamcoverage: bamcoverage
	@$(call TEST_DESC,Memcheck bamcoverage)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMCOVERAGE)/memcheck, \
		rm -rf bamstats-coverage covtmp && \
		$(GRIND) $(REPO_ROOT)/bamcoverage $(TEST_ROOT)/bamstats/400ecoli.bam)
