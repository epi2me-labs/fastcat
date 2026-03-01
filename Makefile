OS := $(shell uname)
ifeq ($(OS), Darwin)
    # mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib -L$(shell brew --prefix curl)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
    ARGP_INC ?= -I$(shell brew --prefix argp-standalone)/include
    CFLAGS ?= -fpic -O3 ${ARGP_INC}
    ZCAT = "gzcat"
    # Try to find clang-tidy and scan-build from brew's LLVM
    LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
    ifneq ($(LLVM_PREFIX),)
        CLANG_TIDY ?= $(LLVM_PREFIX)/bin/clang-tidy
        SCAN_BUILD ?= $(LLVM_PREFIX)/bin/scan-build
    else
        CLANG_TIDY ?= clang-tidy
        SCAN_BUILD ?= scan-build
    endif
else
    ARGP ?=
    ARGP_INC ?=
    CFLAGS ?= -fpic -msse3 -O3 ${ARGP_INC}
    ZCAT = "zcat"
    CLANG_TIDY ?= clang-tidy
    SCAN_BUILD ?= scan-build
endif
CFLAGS += -MMD -MP

VALGRIND ?= valgrind
CLANG_FORMAT ?= clang-format
CLANG_FORMAT_STYLE ?= {BasedOnStyle: WebKit, BreakBeforeBraces: Custom, BraceWrapping: { AfterFunction: true, AfterControlStatement: Never }, PointerAlignment: Right, SortIncludes: true, IncludeBlocks: Regroup, IncludeCategories: [ { Regex: "^<.*>$$", Priority: 1 }, { Regex: "^\\\"(htslib|zlib-ng)/.*\\\"$$", Priority: 2 }, { Regex: "^\\\".*\\\"$$", Priority: 3 } ], ReflowComments: true, AllowShortIfStatementsOnASingleLine: Never, InsertBraces: true, UseTab: Never, IndentWidth: 4, ContinuationIndentWidth: 4, AlignAfterOpenBracket: DontAlign, Cpp11BracedListStyle: false}
DEFLATEVER=1.24

CC ?= gcc
AR ?= ar
STATIC_HTSLIB ?= htslib/libhts.a
HTSLIB_PATCH ?= patches/htslib-aux-parse.patch
HTSLIB_PATCH_STAMP ?= .fastcat_patch_applied
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
    CFLAGS += -DNOTHREADS
endif
FORTIFY_CFLAG ?= -D_FORTIFY_SOURCE=2

SANITIZE ?=
ifeq ($(SANITIZE),ubsan)
    CFLAGS := $(filter-out -O%,$(CFLAGS))
    CFLAGS += -O1 -g -fno-omit-frame-pointer
ifeq ($(OS),Darwin)
    SAN_FLAGS := -fsanitize=address,undefined
else
    SAN_FLAGS := -fsanitize=address,undefined,leak
endif
    EXTRA_CFLAGS += $(SAN_FLAGS) -fno-sanitize-recover=undefined
    EXTRA_LDFLAGS += $(SAN_FLAGS)
    FORTIFY_CFLAG :=
ifeq ($(OS),Darwin)
    # Apple's ASan runtime does not support LeakSanitizer.
    ASAN_OPTIONS ?= halt_on_error=1
else
    ASAN_OPTIONS ?= detect_leaks=1:halt_on_error=1
endif
    UBSAN_OPTIONS ?= print_stacktrace=1:halt_on_error=1
    export ASAN_OPTIONS
    export UBSAN_OPTIONS
else
    # disable asserts in non-sanitizer builds
    CFLAGS += -DNDEBUG
endif

PROFILE ?= 0
PROFILER_LIBDIR ?=
ifeq ($(PROFILE),1)
    EXTRA_CFLAGS += -g3 -fno-omit-frame-pointer
    EXTRA_LIBS += -Wl,--no-as-needed -lprofiler -Wl,--as-needed
    ifneq ($(PROFILER_LIBDIR),)
        EXTRA_LDFLAGS += -L$(PROFILER_LIBDIR) -Wl,-rpath,$(PROFILER_LIBDIR)
    endif
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
# Sanitizer-instrumented binaries are incompatible with valgrind preloading.
ifneq ($(SANITIZE),)
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
RUN_FASTCAT = $(PEPPER) $(REPO_ROOT)/fastcat fastq
RUN_BAMSTATS = $(PEPPER) $(REPO_ROOT)/fastcat bam


.PHONY:
default: fastcat

.PHONY:
test: test_unit test_fastcat test_bamstats test_meta test_bamindex test_fastlint

.PHONY:
test_memory: mem_check_fastcat mem_check_bamstats mem_check_bamindex-build mem_check_bamindex-dump mem_check_bamindex-fetch mem_check_fastlint

.PHONY: format format_check format_staged install_hooks
format:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@find src -type f \( -name '*.c' -o -name '*.h' \) -print0 | xargs -0 $(CLANG_FORMAT) -i \
		-style='$(CLANG_FORMAT_STYLE)'

format_check:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@find src -type f \( -name '*.c' -o -name '*.h' \) -print0 | xargs -0 $(CLANG_FORMAT) --dry-run -Werror \
		-style='$(CLANG_FORMAT_STYLE)'

format_staged:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@tmp=$$(mktemp); \
	trap 'rm -f "$$tmp"' EXIT; \
	git diff --cached --name-only --diff-filter=ACMR -z -- ':(glob)src/**/*.c' ':(glob)src/**/*.h' > "$$tmp"; \
	if [ ! -s "$$tmp" ]; then \
		exit 0; \
	fi; \
	xargs -0 $(CLANG_FORMAT) -i -style='$(CLANG_FORMAT_STYLE)' < "$$tmp"; \
	xargs -0 git add -- < "$$tmp"

.PHONY:
clang_tidy: $(HTSLIB_PATCH_STAMP)
	@command -v $(CLANG_TIDY) >/dev/null 2>&1 || { echo "$(CLANG_TIDY) not found. Install with: apt install clang-tidy"; exit 1; }
	@echo "Running clang-tidy static analysis..."
	@failed_files=$$(mktemp); \
	trap 'rm -f "$$failed_files"' EXIT; \
		find src -type f -name '*.c' ! -path 'src/sdust/*' -print0 | \
		xargs -0 -I{} sh -c '$(CLANG_TIDY) --quiet --warnings-as-errors="*" --header-filter="(^|.*/)src/.*" "$$1" -- -Isrc -Ihtslib -Izlib-ng $(CFLAGS) || echo "$$1" >> "$$2"' sh {} "$$failed_files"; \
	xargs_status=$$?; \
	if [ $$xargs_status -ne 0 ]; then \
		exit $$xargs_status; \
	fi; \
	if [ -s "$$failed_files" ]; then \
		echo "clang-tidy failed for $$(wc -l < "$$failed_files") file(s):"; \
		cat "$$failed_files"; \
		exit 1; \
	fi

.PHONY:
scan_build: clean
	@command -v $(SCAN_BUILD) >/dev/null 2>&1 || { echo "$(SCAN_BUILD) not found. Install with: apt install clang-tools"; exit 1; }
	@echo "Running clang static analyzer..."
	@rm -rf scan-build-results
	@$(SCAN_BUILD) -o scan-build-results \
		--status-bugs \
		--exclude ./htslib \
		--exclude $(CURDIR)/htslib \
		--exclude htslib/ \
		--exclude ./zlib-ng \
		--exclude $(CURDIR)/zlib-ng \
		--exclude zlib-ng/ \
		make fastcat
	@echo "Results saved to scan-build-results/"

.PHONY:
static_analysis: format_check clang_tidy
	@echo "Static analysis checks passed!"

.PHONY:
install_hooks:
	@mkdir -p .githooks
	@chmod +x .githooks/pre-commit
	@git config --local core.hooksPath .githooks
	@echo "Installed git hooks path: .githooks"

.PHONY:
clean:
	rm -rf fastcat bamstats bamindex fastlint bamcoverage libfastcatcore.a $(HTSLIB_PATCH_STAMP) \
		src/fastcat/*.o src/fastcat/*.d \
		src/fastlint/*.o src/fastlint/*.d \
		src/bamstats/*.o src/bamstats/*.d \
		src/bamindex/*.o src/bamindex/*.d \
		src/bamcoverage/*.o src/bamcoverage/*.d \
		src/sdust/*.o src/sdust/*.d \
		src/*.o src/*.d src/sink/*.o src/sink/*.d \
		test/unit/*.o test/unit/*.d test/unity/*.o test/unity/*.d \
		test/unit/unit_tests

.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean

# core library
CORE_OBJS = src/common.o src/fc_log.o src/stats.o src/fc_stats_engine.o src/fastqcomments.o src/fc_read.o src/fc_filter.o src/fc_bam_utils.o src/fc_dust.o src/fc_iter.o src/kseq_gz.o src/fc_writer.o src/fc_catstats_args.o src/fc_program.o src/fc_failures.o src/sink/common.o src/sink/output.o src/sink/perread.o src/sink/perfile.o src/sink/counter.o src/sink/hist.o src/sink/coverage.o src/sink/flagstats.o src/bamcoverage/coverage.o src/sdust/sdust.o src/sdust/kalloc.o src/kh_counter.o src/regiter.o src/fc_input_list.o
CORE_LIB = libfastcatcore.a
FASTCAT_SUBCOMMAND_OBJS = src/fastlint/main.o src/fastlint/args.o src/bamindex/main.o src/bamindex/build_main.o src/bamindex/fetch_main.o src/bamindex/dump_main.o src/bamindex/index.o

# unit tests
UNIT_TEST_SRCS = test/unit/main.c test/unit/test_common.c test/unit/test_stats.c test/unit/test_stats_engine.c test/unit/test_fastqcomments.c test/unit/test_filters.c test/unit/test_rg_parse.c test/unit/test_regiter.c test/unit/test_quality.c test/unit/test_bam_utils.c test/unit/test_qualcalcs.c test/unit/test_iter.c test/unit/test_iter_collect.c test/unit/test_bamcoverage.c test/unity/unity.c
UNITY_DEFS = -DUNITY_INCLUDE_DOUBLE -DUNITY_DOUBLE_PRECISION=1e-6
UNIT_TEST_OBJS = $(UNIT_TEST_SRCS:.c=.o)

test/unit/%.o: test/unit/%.c
	$(CC) -Isrc -Ihtslib -Izlib-ng -Itest/unity $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@
test/unit/%.o: | $(HTSLIB_PATCH_STAMP)

test/unity/%.o: test/unity/%.c
	$(CC) -Itest/unity $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@

test/unit/unit_tests: $(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib -Izlib-ng $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

.PHONY: test_unit
test_unit: test/unit/unit_tests
	@$(call TEST_DESC,Unit tests)
	$(PEPPER) ./test/unit/unit_tests
test_unit: | $(HTSLIB_PATCH_STAMP)


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
		&& make -j 4 libhts.a


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
	$(CC) -Isrc -Ihtslib -Izlib-ng -c -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $< -o $@

src/fastqcomments.o: | $(HTSLIB_PATCH_STAMP)
src/fc_read.o: | $(HTSLIB_PATCH_STAMP)

$(CORE_LIB): $(CORE_OBJS) Makefile
	rm -f $@
	$(AR) rcs $@ $(CORE_OBJS)

-include $(wildcard src/*.d)
-include $(wildcard src/sink/*.d)
-include $(wildcard src/fastcat/*.d)
-include $(wildcard src/fastlint/*.d)
-include $(wildcard src/bamstats/*.d)
-include $(wildcard src/bamindex/*.d)
-include $(wildcard src/bamcoverage/*.d)
-include $(wildcard src/sdust/*.d)
-include $(wildcard test/unit/*.d)
-include $(wildcard test/unity/*.d)

fastcat: src/version.o src/fastcat/main.o $(FASTCAT_SUBCOMMAND_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) zlib-ng/libz.a
	$(CC) -Isrc -Izlib-ng $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

###
# fastcat tests

.PHONY:
test_fastcat: mem_check_fastcat mem_check_fastcat_demultiplex mem_check_fastcat_bam mem_check_fastcat_demultiplex_bam test_fastcat_bam_equivalent test_fastcat_reheader_roundtrip test_fastcat_filters test_fastcat_demultiplex_outputs test_fastcat_demultiplex_split test_fastcat_demux_single_equiv test_fastcat_regression_results test_fastcat_perfile_all_filtered test_fastcat_stdin_data_equiv test_fastcat_fofn_equiv test_fastcat_stdin_requires_dash test_fastcat_bam_input_equiv test_fastcat_bam_input_fofn_equiv test_fastcat_bam_input_stdin_equiv test_fastcat_bam_input_rejects_mapped

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
		$(GRIND) $(REPO_ROOT)/fastcat fastq $(TEST_ROOT)/data/*.fastq.gz > /dev/null)

.PHONY: mem_check_fastcat_bam
mem_check_fastcat_bam: fastcat
	@$(call TEST_DESC,Memcheck fastcat (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/bam, \
		$(GRIND) $(REPO_ROOT)/fastcat fastq $(TEST_ROOT)/data/*.fastq.gz -B > /dev/null)

.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: fastcat
	@$(call TEST_DESC,Memcheck fastcat demultiplex (fastq))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/demux, \
		$(GRIND) $(REPO_ROOT)/fastcat fastq $(TEST_ROOT)/data/*.fastq.gz --demultiplex --output demultiplex > /dev/null)

.PHONY: mem_check_fastcat_demultiplex_bam
mem_check_fastcat_demultiplex_bam: fastcat
	@$(call TEST_DESC,Memcheck fastcat demultiplex (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/memcheck/demux_bam, \
		$(GRIND) $(REPO_ROOT)/fastcat fastq $(TEST_ROOT)/data/*.fastq.gz --demultiplex --output demultiplex -B > /dev/null)

.PHONY: test_fastcat_bam_equivalent
test_fastcat_bam_equivalent: fastcat samtools
	@$(call TEST_DESC,Fastcat BAM equivalence)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/*.fastq.gz --output out_reheader -B | $(REPO_ROOT)/samtools view | $(TEST_ROOT)/sort-sam.py > reheader.sam && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/*.fastq.gz --output out_bam -B | $(REPO_ROOT)/samtools view | $(TEST_ROOT)/sort-sam.py > bam.sam && \
		diff reheader.sam bam.sam)

.PHONY: test_fastcat_reheader_roundtrip
test_fastcat_reheader_roundtrip: fastcat samtools
	@$(call TEST_DESC,Fastcat reheader roundtrip)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/reheader, \
		$(RUN_FASTCAT) $(TEST_ROOT)/reheader/kv.fastq --output out_hist > out.fastq && \
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
		$(RUN_FASTCAT) ../data/*.fastq.gz \
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
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_seen"]!=2||$$h["r_record_ok"]!=0||$$h["r_too_short"]!=2){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)

.PHONY: test_fastcat_stdin_data_equiv
test_fastcat_stdin_data_equiv: fastcat
	@$(call TEST_DESC,Fastcat stdin data matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/stdin_data_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output out_file > out_file.fastq && \
		cat $(TEST_ROOT)/data/bc1.fastq.gz | $(RUN_FASTCAT) - --output out_stdin > out_stdin.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_stdin.fastq,fc_stdin_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,fc_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,fc_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,fc_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)

.PHONY: test_fastcat_fofn_equiv
test_fastcat_fofn_equiv: fastcat
	@$(call TEST_DESC,Fastcat --fofn matches explicit file arguments)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/fofn_equiv, \
		printf "%s\n%s\n" $(TEST_ROOT)/data/bc1.fastq.gz $(TEST_ROOT)/data/bc2.fastq.gz > paths.txt && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz $(TEST_ROOT)/data/bc2.fastq.gz --output out_file > out_file.fastq && \
		cat paths.txt | $(RUN_FASTCAT) --fofn - --output out_fofn > out_fofn.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_fofn.fastq,fc_fofn_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_fofn/runids.tsv,fc_fofn_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_fofn/basecallers.tsv,fc_fofn_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_fofn/per_file.tsv,fc_fofn_perfile) && \
		diff -r out_file/histograms out_fofn/histograms)

.PHONY: test_fastcat_stdin_requires_dash
test_fastcat_stdin_requires_dash: fastcat
	@$(call TEST_DESC,Fastcat stdin modes require explicit '-')
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/stdin_requires_dash, \
		set +e; \
		cat $(TEST_ROOT)/data/bc1.fastq.gz | $(RUN_FASTCAT) --output out_missing_dash > /dev/null 2>&1; \
		rc_data=$$?; \
		printf "%s\n%s\n" $(TEST_ROOT)/data/bc1.fastq.gz $(TEST_ROOT)/data/bc2.fastq.gz | \
			$(RUN_FASTCAT) --fofn --output out_missing_fofn_dash > /dev/null 2>&1; \
		rc_fofn=$$?; \
		set -e; \
		test $$rc_data -ne 0 && test $$rc_fofn -ne 0)

.PHONY: test_fastcat_bam_input_equiv
test_fastcat_bam_input_equiv: fastcat
	@$(call TEST_DESC,Fastcat BAM input (unaligned) matches FASTQ input)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_input_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz $(TEST_ROOT)/data/bc2.fastq.gz --output out_fastq > out_fastq.fastq && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc2.fastq.gz --output tmp2 -B > b2.bam && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam b2.bam --output out_bam > out_bam.fastq && \
		cat out_fastq.fastq \
			| sed -E 's/\t[qQ][sS]:f:[^\t[:space:]]+//g' \
			| $(TEST_ROOT)/sort-fastq.py > fc_bam_input_concat.sorted && \
		cat out_bam.fastq \
			| sed -E 's/\t[qQ][sS]:f:[^\t[:space:]]+//g' \
			| $(TEST_ROOT)/sort-fastq.py > fc_bam_input_concat.expected.sorted && \
		diff fc_bam_input_concat.sorted fc_bam_input_concat.expected.sorted && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/runids.tsv,out_bam/runids.tsv,fc_bam_input_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/basecallers.tsv,out_bam/basecallers.tsv,fc_bam_input_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/per_file.tsv,out_bam/per_file.tsv,fc_bam_input_perfile) && \
		diff -r out_fastq/histograms out_bam/histograms)

.PHONY: test_fastcat_bam_input_fofn_equiv
test_fastcat_bam_input_fofn_equiv: fastcat
	@$(call TEST_DESC,Fastcat BAM input via --fofn matches explicit BAM args)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_input_fofn_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc2.fastq.gz --output tmp2 -B > b2.bam && \
		printf "%s\n%s\n" b1.bam b2.bam > bam_paths.txt && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam b2.bam --output out_file > out_file.fastq && \
		$(RUN_FASTCAT) --input-fmt bam --fofn bam_paths.txt --output out_fofn > out_fofn.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_fofn.fastq,fc_bam_fofn_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_fofn/runids.tsv,fc_bam_fofn_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_fofn/basecallers.tsv,fc_bam_fofn_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_fofn/per_file.tsv,fc_bam_fofn_perfile) && \
		diff -r out_file/histograms out_fofn/histograms)

.PHONY: test_fastcat_bam_input_stdin_equiv
test_fastcat_bam_input_stdin_equiv: fastcat
	@$(call TEST_DESC,Fastcat BAM input from stdin matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_input_stdin_equiv, \
		$(RUN_FASTCAT) $(TEST_ROOT)/data/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam --output out_file > out_file.fastq && \
		cat b1.bam | $(RUN_FASTCAT) --input-fmt bam - --output out_stdin > out_stdin.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_stdin.fastq,fc_bam_stdin_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,fc_bam_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,fc_bam_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,fc_bam_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)

.PHONY: test_fastcat_bam_input_rejects_mapped
test_fastcat_bam_input_rejects_mapped: fastcat
	@$(call TEST_DESC,Fastcat BAM input rejects mapped BAM)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTCAT)/bam_input_rejects_mapped, \
		set +e; \
		$(RUN_FASTCAT) --input-fmt bam $(TEST_ROOT)/bamstats/400ecoli.bam --output out_bad > out.log 2>&1; \
		rc=$$?; \
		set -e; \
		test $$rc -ne 0 && grep -q "only accepts unaligned BAM input" out.log)

###
# bamstats tests

.PHONY: 
test_bamstats: test_bamstats_NM test_bamstats_no_nan_metrics test_bamstats_polya test_bamstats_filters test_bamstats_vs_fastcat_unaligned test_bamstats_vs_fastcat_unaligned_demux test_bamstats_demux_single_equiv test_bamstats_perfile_all_filtered test_bamstats_sample_column mem_check_bamstats
.PHONY: test_bamstats_coverage_demux_layout
test_bamstats: test_bamstats_coverage_demux_layout
.PHONY: test_bamstats_stdin_equiv test_bamstats_stdin_bed_equiv
test_bamstats: test_bamstats_stdin_equiv test_bamstats_stdin_bed_equiv

.PHONY: test_bamstats_demux_single_equiv
test_bamstats_demux_single_equiv: fastcat
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

.PHONY: test_bamstats_sample_column
test_bamstats_sample_column: fastcat
	@$(call TEST_DESC,Bamstats sample_name column only when --sample set)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/sample_column, \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats/400ecoli.bam --output out_nosample > /dev/null && \
		! head -n 1 out_nosample/flagstats.tsv | grep -q "sample_name" && \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats/400ecoli.bam --sample "" --output out_empty > /dev/null && \
		head -n 1 out_empty/flagstats.tsv | grep -q "sample_name" && \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats/400ecoli.bam --sample foo --output out_sample > /dev/null && \
		head -n 1 out_sample/flagstats.tsv | grep -q "sample_name")

.PHONY: test_bamstats_coverage_demux_layout
test_bamstats_coverage_demux_layout: fastcat samtools
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

.PHONY: test_bamstats_stdin_equiv
test_bamstats_stdin_equiv: fastcat samtools
	@$(call TEST_DESC,Bamstats stdin BAM matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/stdin_equiv, \
		cp $(TEST_ROOT)/bamstats/coverage_demux.sam cov.sam && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --output out_file > /dev/null && \
		cat cov.bam | $(RUN_BAMSTATS) - --output out_stdin > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,bs_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,bs_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,bs_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)

.PHONY: test_bamstats_stdin_bed_equiv
test_bamstats_stdin_bed_equiv: fastcat samtools
	@$(call TEST_DESC,Bamstats stdin with BED filter matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/stdin_bed_equiv, \
		cp $(TEST_ROOT)/bamstats/coverage_demux.sam cov.sam && \
		printf "chr1\t0\t75\n" > regions.bed && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --bed regions.bed --output out_file > /dev/null && \
		cat cov.bam | $(RUN_BAMSTATS) - --bed regions.bed --output out_stdin > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,bs_stdin_bed_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,bs_stdin_bed_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,bs_stdin_bed_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)

.PHONY: test_bamstats_NM
test_bamstats_NM: fastcat
	@$(call TEST_DESC,Bamstats NM validation)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/nm, \
		set +e; \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_zeroNM/implausible_nm.sam --output out_bad_implausible > /dev/null 2> err_implausible.log; \
		rc_implausible=$$?; \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_zeroNM/inconsistent_nm.sam --output out_bad_inconsistent > /dev/null 2> err_inconsistent.log; \
		rc_inconsistent=$$?; \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_zeroNM/test.sam --output out_zero > /dev/null 2> err_zero.log; \
		rc_zero=$$?; \
		set -e; \
		test $$rc_implausible -ne 0 && grep -Eq "appears to contain implausible alignment information|zero aligned match/length" err_implausible.log && \
		test $$rc_inconsistent -ne 0 && grep -q "inconsistent NM tag" err_inconsistent.log && \
		test $$rc_zero -eq 0)

.PHONY: test_bamstats_no_nan_metrics
test_bamstats_no_nan_metrics: fastcat
	@$(call TEST_DESC,Bamstats rejects alignments that produce non-finite metrics)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/nan_guard, \
		set +e; $(RUN_BAMSTATS) $(TEST_ROOT)/bamstats_zeroNM/nan_bad.sam --output out_bad_nan > /dev/null 2> err_nan.log; rc=$$?; set -e; \
		test $$rc -ne 0 && grep -q "zero aligned match/length" err_nan.log)

.PHONY: test_bamstats_polya
test_bamstats_polya: fastcat
	@$(call TEST_DESC,Bamstats polya)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/polya, \
		$(RUN_BAMSTATS) $(TEST_ROOT)/bamstats/RCS-100A.bam --poly-a --output out_polya > /dev/null && \
		diff out_polya/histograms/polya.hist $(TEST_ROOT)/bamstats/RCS-100A.bam.polya.hist)

.PHONY: test_bamstats_filters
test_bamstats_filters: fastcat
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
test_bamstats_vs_fastcat_unaligned: fastcat
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
test_bamstats_vs_fastcat_unaligned_demux: fastcat
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
test_bamstats_perfile_all_filtered: fastcat
	@$(call TEST_DESC,Bamstats per-file stats when all reads filtered)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/perfile_filtered, \
		printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n" > short.sam && \
		$(RUN_BAMSTATS) short.sam --tag-name ZZ --tag-value 1 --output out > /dev/null && \
		awk 'NR==1{for(i=1;i<=NF;i++)h[$$i]=i;next} NR==2{ \
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_seen"]!=2||$$h["r_record_ok"]!=0||$$h["r_too_short"]!=0){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)

.PHONY:
mem_check_bamstats: mem_check_bamstats_good mem_check_bamstats_bad mem_check_bamstats_qcfail mem_check_bamstats_duplex mem_check_bamstats_coverage

.PHONY: mem_check_bamstats_good
mem_check_bamstats_good: fastcat
	@$(call TEST_DESC,Memcheck bamstats with good data)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/good, \
		$(GRIND) $(REPO_ROOT)/fastcat bam $(TEST_ROOT)/parse_rg/dna_r10.4.1_e8.2_400bps_hac@v4.3.0.bam > /dev/null)

.PHONY: mem_check_bamstats_bad
mem_check_bamstats_bad: fastcat
	@$(call TEST_DESC,Memcheck bamstats with bad data)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/bad, \
		$(GRIND) $(REPO_ROOT)/fastcat bam $(TEST_ROOT)/parse_rg/bad-ones.bam > /dev/null)

.PHONY: mem_check_bamstats_qcfail
mem_check_bamstats_qcfail: fastcat
	@$(call TEST_DESC,Memcheck bamstats qcfail)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/qcfail, \
		$(GRIND) $(REPO_ROOT)/fastcat bam $(TEST_ROOT)/bamstats/400ecoli-with-qcfail.bam > /dev/null)

.PHONY: mem_check_bamstats_duplex
mem_check_bamstats_duplex: fastcat
	@$(call TEST_DESC,Memcheck bamstats duplex)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/duplex, \
		$(GRIND) $(REPO_ROOT)/fastcat bam $(TEST_ROOT)/bamstats/310dx.bam > /dev/null)

.PHONY: mem_check_bamstats_coverage
mem_check_bamstats_coverage: fastcat
	@$(call TEST_DESC,Memcheck bamstats coverage)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMSTATS)/memcheck/coverage, \
		$(GRIND) $(REPO_ROOT)/fastcat bam --coverage $(TEST_ROOT)/bamstats/400ecoli.bam > /dev/null)

.PHONY:
bamstats_overwrite: fastcat
	@$(call TEST_DESC,Bamstats basecaller file exists)
	@$(PEPPER) cat test/bamstats/400ecoli.bam \
		| ./fastcat bam -u --output bamstats_results - \
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
test_meta_bamstats: fastcat
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
mem_check_bamindex-build: fastcat
	@$(call TEST_DESC,Memcheck bamindex build)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/build, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/fastcat index build 400.bam > /dev/null)

.PHONY: mem_check_bamindex-dump
mem_check_bamindex-dump: fastcat
	@$(call TEST_DESC,Memcheck bamindex dump)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/dump, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/fastcat index build 400.bam > /dev/null && \
		$(GRIND) $(REPO_ROOT)/fastcat index dump 400.bam.bci > /dev/null)

.PHONY: mem_check_bamindex-fetch
mem_check_bamindex-fetch: fastcat
	@$(call TEST_DESC,Memcheck bamindex fetch)
	$(call RUN_IN_TMP,$(TEST_TMP_BAMINDEX)/fetch, \
		cp $(TEST_ROOT)/bamindex/400.bam . && \
		$(GRIND) $(REPO_ROOT)/fastcat index build 400.bam > /dev/null && \
		$(GRIND) $(REPO_ROOT)/fastcat index fetch 400.bam --chunk 5 > /dev/null)

###
# fastlint tests

test_fastlint: mem_check_fastlint

.PHONY:
mem_check_fastlint: fastcat
	@$(call TEST_DESC,Memcheck fastlint)
	$(call RUN_IN_TMP,$(TEST_TMP_FASTLINT)/memcheck, \
		$(GRIND) $(REPO_ROOT)/fastcat lint $(TEST_ROOT)/data/*.fastq.gz > /dev/null)
