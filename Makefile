.DEFAULT_GOAL := default

.PHONY: help
help:  ## Show this help message
	@printf "\n"
	@printf "\033[1mBuild Modes:\033[0m (set with BUILD_MODE=<mode>, default: release)\n"
	@printf "  \033[36mrelease\033[0m      Optimized build for production (-O3, no asserts)\n"
	@printf "  \033[36mdebug\033[0m        Debug build with symbols and asserts (-O0 -g3)\n"
	@printf "  \033[36msanitize\033[0m     Address/undefined/leak sanitizers enabled (-O1)\n"
	@printf "  \033[36mprofile\033[0m      Profiling with gperftools (requires libprofiler)\n"
	@printf "  \033[36mvalgrind\033[0m     Standard build, tests run under valgrind\n"
	@printf "\n"
	@printf "\033[1mBuild Options:\033[0m\n"
	@printf "  \033[36mUSE_DEFLATE=1\033[0m    Use libdeflate compression library\n"
	@printf "  \033[36mNOTHREADS=1\033[0m      Disable threading support\n"
	@printf "\n"
	@printf "\033[1mAvailable Targets:\033[0m\n"
	@grep -hE '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'
	@printf "\n"


# Build directory structure
REPO_ROOT := $(CURDIR)
BUILD_DIR := $(REPO_ROOT)/build
OBJ_DIR := $(BUILD_DIR)/obj
DEP_DIR := $(BUILD_DIR)/dep
FASTCAT_BIN := $(BUILD_DIR)/fastcat
REPORT_DIR := $(BUILD_DIR)/report
REPORT_VENV := $(BUILD_DIR)/report-venv
REPORT_PY := $(REPORT_VENV)/bin/python
REPORT_PIP := $(REPORT_VENV)/bin/pip
REPORT_REQS := tools/requirements.txt
REPORT_MPLCONFIGDIR := $(REPORT_DIR)/mpl-cache
REPORT_EXAMPLE_N ?= 2000
REPORT_DEMO_N ?= 100000
REPORT_DEMO_DURATION_S ?= 172800
DEV_TOOLS_SCRIPT := $(REPO_ROOT)/tools/dev-tools
DEV_TOOLS := python3 $(DEV_TOOLS_SCRIPT)
SUMMARY_JSON_SCHEMA := $(REPO_ROOT)/test/fixtures/misc/summary_json.schema.json
SHELL := /bin/bash
BUILD_MAKEFILES := Makefile

.PHONY: default
default: $(FASTCAT_BIN) ## Build the fastcat binary (default target)

###
# User-customizable variables
# These can be set on the command line or in the environment:
#   EXTRA_CFLAGS, EXTRA_LDFLAGS, EXTRA_LIBS - additional compiler/linker flags
#   HTS_CONF_ARGS - arguments passed to htslib ./configure
#   PROFILER_LIBDIR - custom profiler library directory (with BUILD_MODE=profile)
# Platform-specific defaults are set below

# tool executables
CC ?= gcc
AR ?= ar
VALGRIND ?= valgrind
CLANG_TIDY ?= clang-tidy
SCAN_BUILD ?= scan-build
ZCAT ?= zcat
CLANG_FORMAT ?= clang-format
CLANG_FORMAT_STYLE ?= {BasedOnStyle: WebKit, BreakBeforeBraces: Custom, BraceWrapping: { AfterFunction: true, AfterControlStatement: Never }, PointerAlignment: Right, SortIncludes: true, IncludeBlocks: Regroup, IncludeCategories: [ { Regex: "^<.*>$$", Priority: 1 }, { Regex: "^\\\"(third_party|htslib|zlib-ng)/.*\\\"$$", Priority: 2 }, { Regex: "^\\\".*\\\"$$", Priority: 3 } ], ReflowComments: true, AllowShortIfStatementsOnASingleLine: Never, UseTab: Never, IndentWidth: 4, ContinuationIndentWidth: 4, AlignAfterOpenBracket: DontAlign, AlignEscapedNewlines: DontAlign, Cpp11BracedListStyle: false}

# third-party library versions and paths
DEFLATEVER = 1.24
STATIC_HTSLIB ?= third_party/htslib/libhts.a
HTSLIB_PATCH ?= third_party/patches/htslib-aux-parse.patch
HTSLIB_PATCH_STAMP ?= .fastcat_patch_applied

# macOS: use homebrew packages
OS := $(shell uname)
ifeq ($(OS), Darwin)
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib -L$(shell brew --prefix curl)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
    ARGP_INC ?= -I$(shell brew --prefix argp-standalone)/include
    ZCAT = gzcat
    # Prefer clang-tidy and scan-build from brew's LLVM if available
    LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
    ifneq ($(LLVM_PREFIX),)
        CLANG_TIDY = $(LLVM_PREFIX)/bin/clang-tidy
        SCAN_BUILD = $(LLVM_PREFIX)/bin/scan-build
    endif
endif

# initialize customizable variables if not set
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
PROFILER_LIBDIR ?=
ARGP ?=
ARGP_INC ?=

# base compiler flags (see also BUILD_MODE below)
# Use ?= so conda/environment can override CFLAGS entirely
CFLAGS ?= -fpic ${ARGP_INC}
CFLAGS += -MMD -MP

# detect compiler type (cache result)
CC_IS_CLANG := $(findstring clang,$(shell $(CC) --version 2>/dev/null))

# compiler warnings
ifeq ($(CC_IS_CLANG),)
    # GCC
    # Note: can't use pedantic everywhere because min/max macros trigger:
    #   "ISO C forbids braced-groups within expressions [-Werror=pedantic]"
    WARNINGS = -Werror -Wall -Wextra -Wno-incompatible-pointer-types
else
    # Clang
    # Note: -Wpedantic fine here though
    WARNINGS = -Werror -Wall -Wextra -Wpedantic -Wno-language-extension-token -Wno-gnu-statement-expression -Wno-incompatible-function-pointer-types
endif

###
# Build configuration options
#
# BUILD_MODE: Controls optimization, debugging, and instrumentation (default: release)
#   - release   : Optimized build (-O3), assertions disabled, suitable for production
#   - debug     : No optimization (-O0), debug symbols (-g3), assertions enabled
#   - sanitize  : Address/undefined/leak sanitizers, reduced optimization (-O1)
#   - profile   : Profiling with gperftools, requires libprofiler installed
#   - valgrind  : Optimized build (-O3), tests run under valgrind
#
# USE_DEFLATE=1: Use libdeflate compression library with htslib
# NOTHREADS=1: Disable threading support

BUILD_MODE ?= release

# validate BUILD_MODE
VALID_BUILD_MODES := release debug sanitize profile valgrind
ifeq ($(filter $(BUILD_MODE),$(VALID_BUILD_MODES)),)
    $(error Invalid BUILD_MODE: $(BUILD_MODE). Valid modes: $(VALID_BUILD_MODES))
endif

# disable threading
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
    CFLAGS += -DNOTHREADS
endif

# configure build based on BUILD_MODE
FORTIFY_CFLAG = -D_FORTIFY_SOURCE=2

ifeq ($(BUILD_MODE),debug)
    # debug build: no optimization, full debug symbols, assertions enabled
    CFLAGS += -O0 -g3 -fno-omit-frame-pointer
    # no -DNDEBUG (assertions enabled), no FORTIFY (requires optimization)
    FORTIFY_CFLAG =

else ifeq ($(BUILD_MODE),sanitize)
    # sanitizer build: address/undefined/leak sanitizers, reduced optimization
    CFLAGS += -O1 -g -fno-omit-frame-pointer
ifeq ($(OS),Darwin)
    SAN_FLAGS := -fsanitize=address,undefined
else
    SAN_FLAGS := -fsanitize=address,undefined,leak
endif
    EXTRA_CFLAGS += $(SAN_FLAGS) -fno-sanitize-recover=undefined
    EXTRA_LDFLAGS += $(SAN_FLAGS)
    FORTIFY_CFLAG =
ifeq ($(OS),Darwin)
    ASAN_OPTIONS ?= halt_on_error=1
else
    ASAN_OPTIONS ?= detect_leaks=1:halt_on_error=1
endif
    UBSAN_OPTIONS ?= print_stacktrace=1:halt_on_error=1
    export ASAN_OPTIONS
    export UBSAN_OPTIONS

else ifeq ($(BUILD_MODE),profile)
    # profiling build: optimized with profiling support
    CFLAGS += -O3 -g3 -fno-omit-frame-pointer -DNDEBUG
    EXTRA_LIBS += -Wl,--no-as-needed -lprofiler -Wl,--as-needed
    ifneq ($(PROFILER_LIBDIR),)
        EXTRA_LDFLAGS += -L$(PROFILER_LIBDIR) -Wl,-rpath,$(PROFILER_LIBDIR)
    endif

else ifeq ($(BUILD_MODE),valgrind)
    # valgrind build: optimized build, tests run under valgrind
    CFLAGS += -O3 -DNDEBUG

else ifeq ($(BUILD_MODE),release)
    # release build: optimized, assertions disabled
    CFLAGS += -O3 -DNDEBUG

endif

# add FORTIFY_CFLAG if applicable
ifneq ($(FORTIFY_CFLAG),)
    CFLAGS += $(FORTIFY_CFLAG)
endif

# USE_DEFLATE configuration
HTS_BUILD_FLAGS =
STATIC_LIBDEFLATE =
ifeq ($(USE_DEFLATE), 1)
    # libdeflate paths
    LIBDEF = $(abspath third_party/libdeflate-$(DEFLATEVER))
    LIBDEF_LIB = $(abspath third_party/libdeflate-$(DEFLATEVER)/build)
    # add to build
    STATIC_HTSLIB += third_party/libdeflate-$(DEFLATEVER)/build/libdeflate.a
    STATIC_LIBDEFLATE = third_party/libdeflate-$(DEFLATEVER)/build/libdeflate.a
    HTS_BUILD_FLAGS = CPPFLAGS="-I$(LIBDEF)" LDFLAGS="-L$(LIBDEF_LIB) -Wl,-rpath,$(LIBDEF_LIB)" LIBS="-ldeflate"
endif

# test configuration (valgrind/memcheck)
PEPPER =
GRIND =

# only set up valgrind if not on macOS and not using sanitizers
ifneq ($(OS),Darwin)
    ifneq ($(BUILD_MODE),sanitize)
        GRIND = $(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s
    endif
endif

# enable valgrind for tests in valgrind mode
ifeq ($(BUILD_MODE),valgrind)
    PEPPER = $(GRIND)
endif

# include paths
INCLUDES := -Isrc -Ithird_party/htslib -Ithird_party/zlib-ng

# core library object files
CORE_SRCS := \
	src/core/fc_common.c \
	src/core/fc_log.c \
	src/core/fc_stats.c \
	src/core/fc_quantile.c \
	src/core/fc_kll.c \
	src/core/fc_channel_summary.c \
	src/core/fc_channel_map.c \
	src/core/fc_channel_map_data.c \
	src/core/fc_metric_qscore_2d.c \
	src/core/fc_length_qscore_2d.c \
	src/core/fc_run_summary.c \
	src/core/fc_speed_qscore_2d.c \
	src/core/fc_time_bin.c \
	src/core/fc_stats_engine.c \
	src/core/fc_fastqcomments.c \
	src/core/fc_json.c \
	src/core/fc_report.c \
	src/core/fc_report_arg_table.c \
	src/core/fc_report_template_data.c \
	src/core/fc_bam_rg_cache.c \
	src/core/fc_read.c \
	src/core/fc_filter.c \
	src/core/fc_bam_utils.c \
	src/core/fc_dust.c \
	src/core/fc_iter.c \
	src/core/kseq_gz.c \
	src/core/fc_writer.c \
	src/core/fc_failures.c \
	src/core/kh_counter.c \
	src/core/fc_regiter.c \
	src/core/fc_input_list.c \
	src/core/fc_coverage.c \
	src/sink/common.c \
	src/sink/output.c \
	src/sink/perread.c \
	src/sink/perfile.c \
	src/sink/counter.c \
	src/sink/run_summary.c \
	src/sink/channel.c \
	src/sink/length_qscore_2d.c \
	src/sink/speed_qscore_2d.c \
	src/sink/time_bin.c \
	src/sink/summary_json.c \
	src/sink/hist.c \
	src/sink/coverage.c \
	src/sink/flagstats.c \
	src/vendor/parson/parson.c \
	src/vendor/sdust/sdust.c \
	src/vendor/sdust/kalloc.c

# CLI object files
CLI_SRCS := \
	src/cli/main.c \
	src/cli/catstats/main.c \
	src/cli/catstats/args.c \
	src/cli/bamindex/main.c \
	src/cli/bamindex/build_main.c \
	src/cli/bamindex/fetch_main.c \
	src/cli/bamindex/dump_main.c \
	src/cli/bamindex/index.c \
	src/cli/fastlint/main.c \
	src/cli/fastlint/args.c

# convert source paths to object paths in build/obj/
CORE_OBJS := $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(CORE_SRCS))
CLI_OBJS := $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(CLI_SRCS))
VERSION_OBJ := $(OBJ_DIR)/version.o
CORE_LIB = $(BUILD_DIR)/libfastcatcore.a

###
# build stages

# create build directories
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) $(OBJ_DIR) $(DEP_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(DEP_DIR):
	mkdir -p $(DEP_DIR)

# a patch to htslib to expose some internal APIs we want to hijack
$(HTSLIB_PATCH_STAMP): $(HTSLIB_PATCH)
	@if [ -f "$(HTSLIB_PATCH)" ] && { ! grep -q "sam_aux_parse(" third_party/htslib/htslib/sam.h || ! grep -q "sam_format_aux(" third_party/htslib/htslib/sam.h; }; then \
		echo "Applying htslib patch $(HTSLIB_PATCH)"; \
		cd third_party/htslib && patch -p1 -N < ../../$(HTSLIB_PATCH); \
	fi
	@touch $(HTSLIB_PATCH_STAMP)

third_party/htslib/libhts.a: $(STATIC_LIBDEFLATE) | $(HTSLIB_PATCH_STAMP)
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd third_party/htslib/ \
		&& autoheader \
		&& autoconf \
		&& autoreconf --install \
		&& $(HTS_BUILD_FLAGS) CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure $(HTS_CONF_ARGS) \
		&& $(MAKE) libhts.a

third_party/libdeflate-$(DEFLATEVER)/build/libdeflate.a:
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	mkdir -p third_party && cd third_party && \
	curl -L -o libdeflate-v${DEFLATEVER}.tar.gz https://github.com/ebiggers/libdeflate/archive/refs/tags/v${DEFLATEVER}.tar.gz \
		&& tar -xzf libdeflate-v${DEFLATEVER}.tar.gz \
		&& cd libdeflate-${DEFLATEVER} \
		&& cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
			-DCMAKE_INSTALL_PREFIX="$(pwd)/install" \
		&& cmake --build build -j \

third_party/zlib-ng/zlib.h:
	@echo Configuring zlib-ng
	cd third_party/zlib-ng/ \
		&& CFLAGS="$(CFLAGS) $(EXTRA_CFLAGS)" ./configure --zlib-compat \

third_party/zlib-ng/libz.a: third_party/zlib-ng/zlib.h
	@echo Compiling $(@F)
	cd third_party/zlib-ng/ \
		&& $(MAKE) libz.a

# just for testing
SAMVER=1.21
samtools:
	mkdir -p third_party && cd third_party && \
	curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2 && \
	tar -xjf samtools-${SAMVER}.tar.bz2 && \
	rm samtools-${SAMVER}.tar.bz2 && \
	cd samtools-${SAMVER} && $(MAKE) && \
	cp samtools ../../$@

# outputs to build/obj/ mirroring src/ structure
$(OBJ_DIR)/%.o: src/%.c third_party/zlib-ng/zlib.h $(BUILD_MAKEFILES) | $(OBJ_DIR) $(DEP_DIR)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@))
	$(CC) $(INCLUDES) -c -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@) $< -o $@

# version object
$(VERSION_OBJ): src/version.c third_party/zlib-ng/zlib.h $(BUILD_MAKEFILES) | $(OBJ_DIR) $(DEP_DIR)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@))
	$(CC) $(INCLUDES) -c -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@) $< -o $@

# dependencies on htslib patch for files that need it
$(OBJ_DIR)/core/fc_fastqcomments.o: | $(HTSLIB_PATCH_STAMP)
$(OBJ_DIR)/core/fc_read.o: | $(HTSLIB_PATCH_STAMP)

$(CORE_LIB): $(CORE_OBJS) $(BUILD_MAKEFILES) | $(BUILD_DIR)
	rm -f $@
	$(AR) rcs $@ $(CORE_OBJS)

# include dependency files
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(CORE_OBJS))
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(CLI_OBJS))
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(VERSION_OBJ))

# oh, so thats what you wanted to do
$(FASTCAT_BIN): $(VERSION_OBJ) $(CLI_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) third_party/zlib-ng/libz.a $(BUILD_MAKEFILES) | $(BUILD_DIR)
	@echo "\x1b[1;32mLinking $@\x1b[0m"
	$(CC) -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$(VERSION_OBJ) $(CLI_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) \
		third_party/zlib-ng/libz.a \
		$(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

.PHONY: clean
clean:  ## clean first-party build artifacts
	rm -rf $(BUILD_DIR) $(HTSLIB_PATCH_STAMP) \
		test/unit/*.o test/unit/*.d test/unity/*.o test/unity/*.d \
		test/unit/unit_tests test/tools/*.o

.PHONY: clean_htslib
clean_htslib:
	cd third_party/htslib && [ -f Makefile ] && $(MAKE) clean || true
	rm -f third_party/htslib/config.h third_party/htslib/config.log third_party/htslib/config.status
	rm -f third_party/htslib/config.mk third_party/htslib/htslib.pc third_party/htslib/htslib_static.mk

.PHONY: clean_zlib
clean_zlib:
	cd third_party/zlib-ng && [ -f Makefile ] && $(MAKE) clean || true
	rm -f third_party/zlib-ng/zlib.h third_party/zlib-ng/Makefile

.PHONY: clean_all
clean_all: clean clean_htslib clean_zlib  ## clean all build artifacts including third-party libraries
	rm -rf third_party/libdeflate-* third_party/samtools-*

.PHONY: clang_format
clang_format:  ## run clang-format on first-party source files
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@find src -type f \( -name '*.c' -o -name '*.h' \) ! -path 'src/vendor/*' -print0 | xargs -0 $(CLANG_FORMAT) -i \
		-style='$(CLANG_FORMAT_STYLE)'

.PHONY: clang_format_staged
clang_format_staged:  ## run clang-format on staged source files
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@tmp=$$(mktemp); \
	trap 'rm -f "$$tmp"' EXIT; \
	git diff --cached --name-only --diff-filter=ACMR -z -- ':(glob)src/**/*.c' ':(glob)src/**/*.h' ':(exclude)src/vendor/**' > "$$tmp"; \
	if [ ! -s "$$tmp" ]; then \
		exit 0; \
	fi; \
	xargs -0 $(CLANG_FORMAT) -i -style='$(CLANG_FORMAT_STYLE)' < "$$tmp"; \
	xargs -0 git add -- < "$$tmp"

.PHONY: clang_format_check
clang_format_check:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@find src -type f \( -name '*.c' -o -name '*.h' \) ! -path 'src/vendor/*' -print0 | xargs -0 $(CLANG_FORMAT) --dry-run -Werror \
		-style='$(CLANG_FORMAT_STYLE)'

.PHONY: clang_tidy
clang_tidy: $(HTSLIB_PATCH_STAMP)  ## run clang-tidy static analysis
	@command -v $(CLANG_TIDY) >/dev/null 2>&1 || { echo "$(CLANG_TIDY) not found. Install with: apt install clang-tidy"; exit 1; }
	@echo "Running clang-tidy static analysis..."
	@failed_files=$$(mktemp); \
	trap 'rm -f "$$failed_files"' EXIT; \
		find src -type f -name '*.c' ! -path 'src/vendor/*' -print0 | \
		xargs -0 -I{} sh -c '$(CLANG_TIDY) --quiet --warnings-as-errors="*" --header-filter="(^|.*/)src/.*" "$$1" -- $(INCLUDES) $(CFLAGS) || echo "$$1" >> "$$2"' sh {} "$$failed_files"; \
	xargs_status=$$?; \
	if [ $$xargs_status -ne 0 ]; then \
		exit $$xargs_status; \
	fi; \
	if [ -s "$$failed_files" ]; then \
		echo "clang-tidy failed for $$(wc -l < "$$failed_files") file(s):"; \
		cat "$$failed_files"; \
		exit 1; \
	fi

.PHONY: scan_build
scan_build: clean  ## run scan-build static analysis
	@command -v $(SCAN_BUILD) >/dev/null 2>&1 || { echo "$(SCAN_BUILD) not found. Install with: apt install clang-tools"; exit 1; }
	@echo "Running clang static analyzer..."
	@rm -rf scan-build-results
	@$(SCAN_BUILD) -o scan-build-results \
		--status-bugs \
		--exclude ./third_party/htslib \
		--exclude $(CURDIR)/third_party/htslib \
		--exclude third_party/htslib/ \
		--exclude ./third_party/zlib-ng \
		--exclude $(CURDIR)/third_party/zlib-ng \
		--exclude third_party/zlib-ng/ \
		--exclude ./src/vendor \
		--exclude $(CURDIR)/src/vendor \
		--exclude src/vendor/ \
		$(MAKE) $(FASTCAT_BIN)
	@echo "Results saved to scan-build-results/"

.PHONY: install_hooks
install_hooks:  ## install git hooks for the current repository
	@mkdir -p .githooks
	@chmod +x .githooks/pre-commit
	@git config --local core.hooksPath .githooks
	@echo "Installed git hooks path: .githooks"


###
# Developer report scaffold
#

$(REPORT_DIR):
	mkdir -p $(REPORT_DIR)

$(REPORT_MPLCONFIGDIR): | $(REPORT_DIR)
	mkdir -p $(REPORT_MPLCONFIGDIR)

$(REPORT_PY):
	python3 -m venv $(REPORT_VENV)

$(REPORT_VENV)/.deps: $(REPORT_REQS) $(REPORT_PY) | $(REPORT_MPLCONFIGDIR)
	$(REPORT_PIP) install -r $(REPORT_REQS)
	@touch $@

.PHONY: report_venv
report_venv: $(REPORT_VENV)/.deps  ## create the Python report venv

$(REPORT_DIR)/template_example.fastq: $(DEV_TOOLS_SCRIPT) tools/workflow_glue/simulate_reads.py | $(REPORT_DIR)
	$(DEV_TOOLS) simulate_reads \
		--n $(REPORT_EXAMPLE_N) \
		--shape 2.6 \
		--scale 30000 \
		--q1 10 \
		--q2 20 \
		--q1-weight 0.55 \
		--seed 42 \
		> $@

$(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv: $(FASTCAT_BIN) $(REPORT_DIR)/template_example.fastq | $(REPORT_DIR)
	rm -rf $(REPORT_DIR)/template_fastcat
	$(FASTCAT_BIN) fastq $(REPORT_DIR)/template_example.fastq \
		--sample __FASTCAT_SAMPLE_ALIAS_HTML__ \
		--output $(REPORT_DIR)/template_fastcat \
		> /dev/null

$(REPORT_DIR)/template_fastcat/histograms/length.hist: $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv
	test -s $@

$(REPORT_DIR)/template_fastcat/time_bin_summary.tsv: $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv
	test -s $@

$(REPORT_DIR)/template_fastcat/speed_qscore_2d.tsv: $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv
	test -s $@

$(REPORT_DIR)/template_fastcat/channel_summary.tsv: $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv
	test -s $@

$(REPORT_DIR)/template_inputs/metadata.json: $(DEV_TOOLS_SCRIPT) tools/workflow_glue/write_wf_report_inputs.py | $(REPORT_DIR)
	$(DEV_TOOLS) write_wf_report_inputs --template --out-dir $(REPORT_DIR)/template_inputs

$(REPORT_DIR)/template_inputs/run_summary.tsv: $(REPORT_DIR)/template_inputs/metadata.json
	test -s $@

$(REPORT_DIR)/example_inputs/metadata.json: $(DEV_TOOLS_SCRIPT) tools/workflow_glue/write_wf_report_inputs.py | $(REPORT_DIR)
	$(DEV_TOOLS) write_wf_report_inputs --out-dir $(REPORT_DIR)/example_inputs

$(REPORT_DIR)/example_inputs/run_summary.tsv: $(REPORT_DIR)/example_inputs/metadata.json
	test -s $@

$(REPORT_DIR)/template_wf_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/template_fastcat/histograms/length.hist $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv $(REPORT_DIR)/template_fastcat/speed_qscore_2d.tsv $(REPORT_DIR)/template_fastcat/time_bin_summary.tsv $(REPORT_DIR)/template_fastcat/channel_summary.tsv $(REPORT_DIR)/template_inputs/metadata.json $(REPORT_DIR)/template_inputs/run_summary.tsv $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--stats $(REPORT_DIR)/template_fastcat/histograms \
		--run-summary $(REPORT_DIR)/template_inputs/run_summary.tsv \
		--length-qscore-2d $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv \
		--speed-qscore-2d $(REPORT_DIR)/template_fastcat/speed_qscore_2d.tsv \
		--time-bin-summary $(REPORT_DIR)/template_fastcat/time_bin_summary.tsv \
		--channel-summary $(REPORT_DIR)/template_fastcat/channel_summary.tsv \
		--metadata $(REPORT_DIR)/template_inputs/metadata.json \
		--versions $(REPORT_DIR)/template_inputs/versions \
		--params $(REPORT_DIR)/template_inputs/params.json \
		--wf_version __FASTCAT_VERSION_HTML__

$(REPORT_DIR)/example_wf_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/template_fastcat/histograms/length.hist $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv $(REPORT_DIR)/template_fastcat/speed_qscore_2d.tsv $(REPORT_DIR)/template_fastcat/time_bin_summary.tsv $(REPORT_DIR)/template_fastcat/channel_summary.tsv $(REPORT_DIR)/example_inputs/metadata.json $(REPORT_DIR)/example_inputs/run_summary.tsv $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--stats $(REPORT_DIR)/template_fastcat/histograms \
		--run-summary $(REPORT_DIR)/example_inputs/run_summary.tsv \
		--length-qscore-2d $(REPORT_DIR)/template_fastcat/length_qscore_2d.tsv \
		--speed-qscore-2d $(REPORT_DIR)/template_fastcat/speed_qscore_2d.tsv \
		--time-bin-summary $(REPORT_DIR)/template_fastcat/time_bin_summary.tsv \
		--channel-summary $(REPORT_DIR)/template_fastcat/channel_summary.tsv \
		--metadata $(REPORT_DIR)/example_inputs/metadata.json \
		--versions $(REPORT_DIR)/example_inputs/versions \
		--params $(REPORT_DIR)/example_inputs/params.json \
		--wf_version dev

$(REPORT_DIR)/time_preview_fastcat/time_bin_summary.tsv: $(FASTCAT_BIN) $(REPORT_DIR)/template_example.fastq | $(REPORT_DIR)
	rm -rf $(REPORT_DIR)/time_preview_fastcat
	$(FASTCAT_BIN) fastq $(REPORT_DIR)/template_example.fastq \
		--sample __FASTCAT_SAMPLE_ALIAS_HTML__ \
		--time-bin-size 60 \
		--output $(REPORT_DIR)/time_preview_fastcat \
		> /dev/null

$(REPORT_DIR)/time_plots_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/time_preview_fastcat/time_bin_summary.tsv $(REPORT_DIR)/time_preview_fastcat/channel_summary.tsv $(REPORT_DIR)/example_inputs/metadata.json $(REPORT_DIR)/example_inputs/run_summary.tsv $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--stats $(REPORT_DIR)/time_preview_fastcat/histograms \
		--run-summary $(REPORT_DIR)/example_inputs/run_summary.tsv \
		--length-qscore-2d $(REPORT_DIR)/time_preview_fastcat/length_qscore_2d.tsv \
		--speed-qscore-2d $(REPORT_DIR)/time_preview_fastcat/speed_qscore_2d.tsv \
		--time-bin-summary $(REPORT_DIR)/time_preview_fastcat/time_bin_summary.tsv \
		--channel-summary $(REPORT_DIR)/time_preview_fastcat/channel_summary.tsv \
		--metadata $(REPORT_DIR)/example_inputs/metadata.json \
		--versions $(REPORT_DIR)/example_inputs/versions \
		--params $(REPORT_DIR)/example_inputs/params.json \
		--wf_version dev

.PHONY: report_time_plots_preview
report_time_plots_preview: $(REPORT_DIR)/time_plots_report.html  ## build Python-only report preview with time-bin plots

.PHONY: report_template
report_template: report_vendor_template  ## build the wf-template-derived HTML report template

.PHONY: report_vendor_template
report_vendor_template: $(REPORT_DIR)/template_wf_report.html $(DEV_TOOLS_SCRIPT) tools/workflow_glue/vendor_report_template.py  ## vendor report template into C assets
	$(DEV_TOOLS) vendor_report_template \
		--input $(REPORT_DIR)/template_wf_report.html \
		--template-output $(REPORT_DIR)/fastcat_report_template.html \
		--header src/core/fc_report_template_data.h \
		--source src/core/fc_report_template_data.c \
		--symbol fc_report_template_html

.PHONY: report_scaffold
report_scaffold: report_vendor_template  ## rebuild the developer HTML report scaffold


###
# Testing!
#

.PHONY: test
test: test_functional test_unit  ## run all tests

# unit tests use unity framework

UNIT_TEST_SRCS = test/unit/main.c $(sort $(wildcard test/unit/test_*.c)) test/unity/unity.c
UNITY_DEFS = -DUNITY_INCLUDE_DOUBLE -DUNITY_DOUBLE_PRECISION=1e-6
UNIT_TEST_OBJS = $(UNIT_TEST_SRCS:.c=.o)
.SECONDARY: $(UNIT_TEST_OBJS) $(CORE_LIB)

test/unit/%.o: test/unit/%.c $(BUILD_MAKEFILES)
	$(CC) $(INCLUDES) -Itest/unity $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@
test/unit/%.o: | $(HTSLIB_PATCH_STAMP)

test/unity/%.o: test/unity/%.c $(BUILD_MAKEFILES)
	$(CC) -Itest/unity $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(UNITY_DEFS) -c $< -o $@

test/unit/unit_tests: $(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) $(BUILD_MAKEFILES)
	$(CC) $(INCLUDES) $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

-include $(UNIT_TEST_SRCS:.c=.d)

.PHONY: test_unit
test_unit: test/unit/unit_tests  ## run unit tests
	@$(call TEST_DESC,Unit tests)
	$(PEPPER) ./test/unit/unit_tests
test_unit: | $(HTSLIB_PATCH_STAMP)


# integration tests
# - Group-to-group registrations are defined here at the top
# - Individual test-to-leaf-group registrations are inline next to each test implementation

# Declare all test group targets as PHONY
.PHONY: test test_functional test_memory test_unit demo demo_ubam test_report_missing_metadata
.PHONY: test_fastcat test_fastcat_memory test_fastcat_equivalence test_fastcat_demux test_fastcat_inputs test_fastcat_filters test_fastcat_output
.PHONY: test_bamstats test_bamstats_memory test_bamstats_core test_bamstats_equivalence test_bamstats_coverage test_bamstats_inputs
.PHONY: test_bamindex test_bamindex_memory
.PHONY: test_fastlint test_fastlint_memory
.PHONY: test_meta

.PHONY: test_functional
test_functional: test_fastcat test_bamstats test_bamindex test_fastlint test_meta demo  ## run integration tests

.PHONY: test_memory
test_memory: test_fastcat_memory test_bamstats_memory test_bamindex_memory test_fastlint_memory

# test_fastcat sub-groups
.PHONY: test_fastcat
test_fastcat: test_fastcat_memory test_fastcat_demux test_fastcat_equivalence test_fastcat_filters test_fastcat_output test_fastcat_inputs

# test_bamstats sub-groups
.PHONY: test_bamstats
test_bamstats: test_bamstats_memory test_bamstats_core test_bamstats_coverage test_bamstats_inputs test_bamstats_equivalence

# test_bamindex sub-groups
.PHONY: test_bamindex
test_bamindex: test_bamindex_memory

# test_fastlint sub-groups
.PHONY: test_fastlint
test_fastlint: test_fastlint_memory

# test_meta is a leaf group

# various helper macros for tests
TEST_ROOT := $(REPO_ROOT)/test
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
COMPARE_TABLES_SKIP_THIRD = \
	head -n 1 $(1) > $(3).header && \
	head -n 1 $(2) > $(3).expected.header && \
	diff $(3).header $(3).expected.header && \
	tail -n +2 $(1) | cut -f1-2,4- | LC_ALL=C sort > $(3).sorted && \
	tail -n +2 $(2) | cut -f1-2,4- | LC_ALL=C sort > $(3).expected.sorted && \
	diff $(3).sorted $(3).expected.sorted
COMPARE_FASTQ = \
	case "$(1)" in *.gz) gunzip -c "$(1)" ;; *) cat "$(1)" ;; esac \
		| $(DEV_TOOLS) sort_fastq > $(3).sorted && \
	case "$(2)" in *.gz) gunzip -c "$(2)" ;; *) cat "$(2)" ;; esac \
		| $(DEV_TOOLS) sort_fastq > $(3).expected.sorted && \
	diff $(3).sorted $(3).expected.sorted
# Validate and compare JSON files, normalizing paths and optionally barcode IDs.
# Args: $(1)=actual $(2)=expected $(3)=optional: mask_barcode
COMPARE_JSON = $(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) --expected $(2) --path-root $(REPO_ROOT) $(if $(filter mask_barcode,$(3)),--mask-barcode) $(1)
TEST_DESC = printf '\n%s\n' "$(1)"
RUN_FASTCAT = $(PEPPER) $(FASTCAT_BIN) fastq
RUN_BAMSTATS = $(PEPPER) $(FASTCAT_BIN) bam

TEST_FIXTURES := $(TEST_ROOT)/fixtures
TEST_FIXTURES_FASTQ := $(TEST_FIXTURES)/fastq
TEST_FIXTURES_BAM := $(TEST_FIXTURES)/bam
TEST_EXPECTED := $(TEST_ROOT)/expected
TEST_EXPECTED_FASTCAT := $(TEST_EXPECTED)/fastcat
TEST_EXPECTED_BAMSTATS := $(TEST_EXPECTED)/bamstats
TEST_EXPECTED_METADATA := $(TEST_EXPECTED)/metadata
TEST_EXPECTED_DEMULTIPLEX := $(TEST_EXPECTED_FASTCAT)/demultiplex
TEST_EXPECTED_SUMMARY_JSON := $(TEST_EXPECTED_FASTCAT)
TEST_TMP_DIR ?= $(TEST_ROOT)/tmp
comma := ,

# Test demux single barcode equivalence
# This is used with both `fastcat fastq` and `fastcat bam` to provide the same test
# $(1) = setup commands to run before comparison
# $(2) = expected files list
# $(3) = test name prefix
# $(4) = per_read.tsv comparison macro (COMPARE_TABLES or COMPARE_TABLES_SKIP_FIRST)
define TEST_DEMUX_SINGLE_EQUIV
	$(1) && \
	$(call COMPARE_TABLES_SKIP_FIRST,plain_out/per_file.tsv,demux_out/per_file.tsv,$(3)_perfile) && \
	$(call COMPARE_TABLES_SKIP_FIRST,plain_out/runids.tsv,demux_out/runids.tsv,$(3)_runids_top) && \
	$(call COMPARE_TABLES_SKIP_FIRST,plain_out/basecallers.tsv,demux_out/basecallers.tsv,$(3)_basecallers_top) && \
	for f in $(2); do \
		echo "Checking $$f..." && \
		if [[ "$$f" == *.hist ]]; then \
			diff plain_out/$$f demux_out/barcode0001/$$f || exit 1; \
		elif [[ "$$f" == "per_read.tsv" ]]; then \
			$(call $(4),plain_out/$$f,demux_out/barcode0001/$$f,$(3)_$$f) || exit 1; \
		elif [[ "$$f" == "run_summary.tsv" ]]; then \
			$(call COMPARE_TABLES_SKIP_THIRD,plain_out/$$f,demux_out/barcode0001/$$f,$(3)_$$f) || exit 1; \
		elif [[ "$$f" == "summary.json" ]]; then \
			$(call COMPARE_JSON,plain_out/$$f,demux_out/barcode0001/$$f,mask_barcode) || exit 1; \
		else \
			$(call COMPARE_TABLES_SKIP_FIRST,plain_out/$$f,demux_out/barcode0001/$$f,$(3)_$$f) || exit 1; \
		fi; \
	done && \
	echo "Checking no unexpected files..." && \
	printf "%s\n" $(2) | LC_ALL=C sort > expected_files.txt && \
	find plain_out -type f ! -name "per_file.tsv" | sed 's|^plain_out/||' | LC_ALL=C sort > plain_files.txt && \
	diff expected_files.txt plain_files.txt && \
	find demux_out/barcode0001 -type f ! -name "*.fastq.gz" | sed 's|^demux_out/barcode0001/||' | LC_ALL=C sort > demux_files.txt && \
	diff expected_files.txt demux_files.txt
endef


###
# meta data tests (both fastcat and bamstats)

# Test metadata extraction (runids and basecallers)
# $(1) = command to run (e.g., $(RUN_FASTCAT) or $(RUN_BAMSTATS))
# $(2) = fixtures path (e.g., $(TEST_FIXTURES_FASTQ)/metadata)
# $(3) = symlink name (e.g., parse_rd or parse_rg)
# $(4) = file glob pattern (e.g., *.fastq or *.bam)
# $(5) = expected files path (e.g., $(TEST_EXPECTED_METADATA))
define TEST_METADATA
	ln -s $(2) $(3) && mkdir -p work && \
	pushd work > /dev/null && \
	set -e; \
	for i in ../$(3)/$(4); do \
		echo $$i; \
		filename=$$(basename $$i); \
		$(1) $$i --output out --per-read > /dev/null; \
		diff out/runids.tsv $(5)/$$filename.runids || exit 1; \
		diff out/basecallers.tsv $(5)/$$filename.callers || exit 1; \
		rm -rf out; \
	done; \
	popd > /dev/null
endef
# Expected files to compare (excluding per_file.tsv which only exists at top level)
FASTCAT_EXPECTED_FILES = runids.tsv basecallers.tsv per_read.tsv run_summary.tsv summary.json channel_summary.tsv length_qscore_2d.tsv speed_qscore_2d.tsv time_bin_summary.tsv histograms/length.hist histograms/quality.hist
BAMSTATS_EXPECTED_FILES = runids.tsv basecallers.tsv per_read.tsv run_summary.tsv summary.json channel_summary.tsv length_qscore_2d.tsv speed_qscore_2d.tsv time_bin_summary.tsv histograms/length.hist histograms/quality.hist histograms/accuracy.hist histograms/coverage.hist histograms/length.unmap.hist histograms/quality.unmap.hist flagstats.tsv


###
# report tests

.PHONY: demo
demo: $(FASTCAT_BIN)  ## generate a synthetic FASTQ run and HTML report demo
	@$(call TEST_DESC,HTML report demo)
	rm -rf $(TEST_TMP_DIR)/report_html
	mkdir -p $(TEST_TMP_DIR)/report_html
	set -o pipefail; $(DEV_TOOLS) simulate_reads \
		--n $(REPORT_DEMO_N) \
		--run-duration-s $(REPORT_DEMO_DURATION_S) \
		--sample-id simulated-data \
		--max-len 1000000 \
		--seed 7 \
		| gzip -c > $(TEST_TMP_DIR)/report_html/reads.fastq.gz
	$(RUN_FASTCAT) $(TEST_TMP_DIR)/report_html/reads.fastq.gz \
		--output $(TEST_TMP_DIR)/report_html/fastcat_out \
		--html-report \
		--json-summary \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/length_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/speed_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/time_bin_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/channel_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/run_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/summary.json
	# TODO: these HTML content greps are brittle; prefer structural checks or remove low-value text assertions
	grep -q "Fastcat Sequencing QC Report" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Quality diagnostics" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Read length" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Translocation speed" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Sequencing over time" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Cumulative base yield" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Read count" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Channel occupancy" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "__FASTCAT_" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "fc_patch_fastcat_report" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Run configuration" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "Software versions" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "Workflow parameters" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "fastcat workflow provided by Oxford Nanopore Technologies" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "nextflow workflow" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Bokeh" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	$(RUN_FASTCAT) $(TEST_TMP_DIR)/report_html/reads.fastq.gz \
		--output $(TEST_TMP_DIR)/report_html/demux_out \
		--demultiplex \
		--html-report \
		--json-summary \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/fastcat_report.html
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/length_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/speed_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/time_bin_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/channel_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/run_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/summary.json
	rm -f $(TEST_TMP_DIR)/report_html/reads.fastq.gz

.PHONY: demo_ubam
demo_ubam: $(FASTCAT_BIN) samtools  ## generate a synthetic unaligned BAM run and HTML report demo
	@$(call TEST_DESC,HTML report uBAM demo)
	rm -rf $(TEST_TMP_DIR)/report_html_ubam
	mkdir -p $(TEST_TMP_DIR)/report_html_ubam
	set -o pipefail; $(DEV_TOOLS) simulate_reads \
		--n $(REPORT_DEMO_N) \
		--run-duration-s $(REPORT_DEMO_DURATION_S) \
		--sample-id simulated-data \
		--max-len 1000000 \
		--seed 7 \
		| gzip -c > $(TEST_TMP_DIR)/report_html_ubam/reads.fastq.gz
	$(REPO_ROOT)/samtools import -T '*' -O BAM -@ 1 \
		-o $(TEST_TMP_DIR)/report_html_ubam/reads.bam \
		$(TEST_TMP_DIR)/report_html_ubam/reads.fastq.gz
	$(RUN_BAMSTATS) $(TEST_TMP_DIR)/report_html_ubam/reads.bam \
		--output $(TEST_TMP_DIR)/report_html_ubam/bam_out \
		--html-report \
		--json-summary \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/length_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/speed_qscore_2d.tsv
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/time_bin_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/channel_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/run_summary.tsv
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/summary.json
	grep -q $$'\\tunaligned bam\\t' $(TEST_TMP_DIR)/report_html_ubam/bam_out/run_summary.tsv
	grep -q "Input type" $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	grep -q "unaligned bam" $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	rm -f $(TEST_TMP_DIR)/report_html_ubam/reads.fastq.gz

.PHONY: test_report_missing_metadata
test_report_missing_metadata: $(FASTCAT_BIN)
	@$(call TEST_DESC,HTML report null-state with missing metadata)
	rm -rf $(TEST_TMP_DIR)/report_missing_metadata
	mkdir -p $(TEST_TMP_DIR)/report_missing_metadata
	set -o pipefail; $(DEV_TOOLS) simulate_reads \
		--n 200 \
		--seed 11 \
		--omit-metadata \
		| gzip -c > $(TEST_TMP_DIR)/report_missing_metadata/reads.fastq.gz
	$(RUN_FASTCAT) $(TEST_TMP_DIR)/report_missing_metadata/reads.fastq.gz \
		--output $(TEST_TMP_DIR)/report_missing_metadata/out \
		--html-report \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No cumulative base yield through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No base yield through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No cumulative read count through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No read count through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No read length through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No quality through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No translocation speed through time data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No translocation speed vs quality data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	grep -q "No channel occupancy data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	! grep -q "No read length vs quality data available" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	! grep -q "__FASTCAT_" $(TEST_TMP_DIR)/report_missing_metadata/out/fastcat_report.html
	rm -f $(TEST_TMP_DIR)/report_missing_metadata/reads.fastq.gz


###
# fastcat tests

.PHONY: test_fastcat_demux_single_equiv
test_fastcat_demux_single_equiv: $(FASTCAT_BIN) $(REPORT_VENV)/.deps
	@$(call TEST_DESC,Fastcat demux vs non-demux (single barcode))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/demux_single, \
		$(call TEST_DEMUX_SINGLE_EQUIV,\
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output plain_out --per-read --json-summary > /dev/null && \
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output demux_out --per-read --demultiplex --json-summary > /dev/null,\
			$(FASTCAT_EXPECTED_FILES),\
			demux_single,\
			COMPARE_TABLES_SKIP_FIRST))
test_fastcat_demux: test_fastcat_demux_single_equiv


.PHONY: test_meta_fastcat
test_meta_fastcat: $(FASTCAT_BIN)
	@$(call TEST_DESC,Metadata fastcat)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/meta/fastq, \
		$(call TEST_METADATA,$(RUN_FASTCAT),$(TEST_FIXTURES_FASTQ)/metadata,parse_rd,*.fastq,$(TEST_EXPECTED_METADATA)))
test_meta: test_meta_fastcat


.PHONY: test_fastcat_bam_equivalent
test_fastcat_bam_equivalent: $(FASTCAT_BIN) samtools
	@$(call TEST_DESC,Fastcat BAM equivalence)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_equiv, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz --output out_reheader -B | $(REPO_ROOT)/samtools view | $(DEV_TOOLS) sort_sam > reheader.sam && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz --output out_bam -B | $(REPO_ROOT)/samtools view | $(DEV_TOOLS) sort_sam > bam.sam && \
		diff reheader.sam bam.sam)
test_fastcat_equivalence: test_fastcat_bam_equivalent


.PHONY: test_fastcat_reheader_roundtrip
test_fastcat_reheader_roundtrip: $(FASTCAT_BIN) samtools
	@$(call TEST_DESC,Fastcat reheader roundtrip)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/reheader, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/reheader/kv.fastq --output out_hist > out.fastq && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/reheader/kv.fastq --output out_orig -B | $(REPO_ROOT)/samtools view > orig.sam && \
		$(RUN_FASTCAT) out.fastq --output out_reheader -B | $(REPO_ROOT)/samtools view > reheader.sam && \
		diff orig.sam reheader.sam)
test_fastcat_equivalence: test_fastcat_reheader_roundtrip


.PHONY: test_fastcat_length_quality_filters
test_fastcat_length_quality_filters: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat filters)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/filters, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/filters/fastcat_filters.fastq -a 5 -b 8 -q 20 --output out_filters > out.fastq && \
		grep -q "^@read_ok" out.fastq && \
		! grep -q "^@read_short" out.fastq && \
		! grep -q "^@read_long" out.fastq && \
		! grep -q "^@read_lowq" out.fastq)
test_fastcat_filters: test_fastcat_length_quality_filters


.PHONY: test_fastcat_demultiplex_outputs
test_fastcat_demultiplex_outputs: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat demultiplex outputs)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/demultiplex/outputs, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc0.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz \
			$(TEST_FIXTURES_FASTQ)/basic/bcEmpty.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bcMixed.fastq.gz --demultiplex --output out_demux > /dev/null && \
		$(call COMPARE_FASTQ,out_demux/unclassified/unclassified.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/unclassified.fastq.gz,demux_unclassified) && \
		$(call COMPARE_FASTQ,out_demux/barcode0001/barcode0001.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001.fastq.gz,demux_barcode0001) && \
		$(call COMPARE_FASTQ,out_demux/barcode0002/barcode0002.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0002.fastq.gz,demux_barcode0002) && \
		test -f out_demux/unclassified/histograms/length.hist && \
		test -f out_demux/barcode0001/histograms/length.hist && \
		test -f out_demux/barcode0001/histograms/quality.hist && \
		test -f out_demux/unclassified/length_qscore_2d.tsv && \
		test -f out_demux/unclassified/speed_qscore_2d.tsv && \
		test -f out_demux/barcode0001/length_qscore_2d.tsv && \
		test -f out_demux/barcode0001/speed_qscore_2d.tsv && \
		test -f out_demux/barcode0002/length_qscore_2d.tsv && \
		test -f out_demux/barcode0002/speed_qscore_2d.tsv && \
		test -f out_demux/unclassified/time_bin_summary.tsv && \
		test -f out_demux/barcode0001/time_bin_summary.tsv && \
		test -f out_demux/barcode0002/time_bin_summary.tsv && \
		test -f out_demux/unclassified/run_summary.tsv && \
		test -f out_demux/barcode0001/run_summary.tsv && \
		test -f out_demux/barcode0002/run_summary.tsv && \
		test -f out_demux/unclassified/channel_summary.tsv && \
		test -f out_demux/barcode0001/channel_summary.tsv && \
		test -f out_demux/barcode0002/channel_summary.tsv)
test_fastcat_demux: test_fastcat_demultiplex_outputs


.PHONY: test_fastcat_demultiplex_split
test_fastcat_demultiplex_split: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat demultiplex split files)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/demultiplex/split, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --demultiplex --output out_split --reads-per-file 4 --json-summary > /dev/null && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0000.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0000.fastq.gz,demux_split_0000) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0001.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0001.fastq.gz,demux_split_0001) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0002.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0002.fastq.gz,demux_split_0002) && \
		test -f out_split/barcode0001/histograms/length.hist && \
		test -f out_split/barcode0001/length_qscore_2d.tsv && \
		test -f out_split/barcode0001/speed_qscore_2d.tsv && \
		test -f out_split/barcode0001/time_bin_summary.tsv && \
		test -f out_split/barcode0001/run_summary.tsv && \
		test -f out_split/barcode0001/channel_summary.tsv && \
		test -f out_split/barcode0001/summary.json)
test_fastcat_demux: test_fastcat_demultiplex_split


.PHONY: test_fastcat_length_qscore_2d
test_fastcat_length_qscore_2d: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat length_qscore_2d output shape and totals)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/length_qscore_2d, \
		printf "@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\nIIII\n@r3\nGGGG\n+\n5555\n" > tiny.fastq && \
		$(RUN_FASTCAT) tiny.fastq --output out > /dev/null && \
		test -f out/length_qscore_2d.tsv && \
		test $$(wc -l < out/length_qscore_2d.tsv) -eq 120001 && \
		awk 'BEGIN { eps = 1e-6; min_qs = 1e9; max_qe = -1e9; } \
			NR == 1 { next; } \
			{ \
				key = $$3 FS $$4; bins[key] = 1; \
				if ($$3 < min_qs) min_qs = $$3; \
				if ($$4 > max_qe) max_qe = $$4; \
				if (sqrt((($$4 - $$3) - 0.5) * (($$4 - $$3) - 0.5)) > eps) exit 1; \
				if ($$5 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/) exit 1; \
				v = $$5 + 0.0; \
				if (v < -eps) exit 1; \
			} \
			END { \
				if (length(bins) != 120) exit 1; \
				if (sqrt(min_qs * min_qs) > eps) exit 1; \
				if (sqrt((max_qe - 60.0) * (max_qe - 60.0)) > eps) exit 1; \
			}' out/length_qscore_2d.tsv)
test_fastcat_output: test_fastcat_length_qscore_2d


.PHONY: test_fastcat_speed_qscore_2d
test_fastcat_speed_qscore_2d: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat speed_qscore_2d output shape and totals)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/speed_qscore_2d, \
		printf "@r1 du=1.0\nACGT\n+\n!!!!\n@r2 du=2.0\nAAAAAA\n+\nIIIIII\n@r3 du=4.0\nGGGGGGGG\n+\n55555555\n@r4\nTTTT\n+\nIIII\n" > tiny.fastq && \
		$(RUN_FASTCAT) tiny.fastq --output out > /dev/null && \
		test -f out/speed_qscore_2d.tsv && \
		test $$(wc -l < out/speed_qscore_2d.tsv) -eq 120001 && \
		awk 'BEGIN { eps = 1e-6; min_qs = 1e9; max_qe = -1e9; } \
			NR == 1 { next; } \
			{ \
				key = $$3 FS $$4; bins[key] = 1; \
				if ($$3 < min_qs) min_qs = $$3; \
				if ($$4 > max_qe) max_qe = $$4; \
				if (sqrt((($$4 - $$3) - 0.5) * (($$4 - $$3) - 0.5)) > eps) exit 1; \
				if ($$5 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/) exit 1; \
				v = $$5 + 0.0; \
				if (v < -eps) exit 1; \
				sum += v; \
			} \
			END { \
				if (length(bins) != 120) exit 1; \
				if (sqrt(min_qs * min_qs) > eps) exit 1; \
				if (sqrt((max_qe - 60.0) * (max_qe - 60.0)) > eps) exit 1; \
				if (sqrt((sum - 3.0) * (sum - 3.0)) > 1e-3) exit 1; \
			}' out/speed_qscore_2d.tsv)
test_fastcat_output: test_fastcat_speed_qscore_2d


.PHONY: test_fastcat_time_bin_summary
test_fastcat_time_bin_summary: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat time_bin_summary output shape and totals)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/time_bin_summary, \
		printf "@r1 start_time=1970-01-01T00:00:00Z du=2.0\nACGT\n+\nIIII\n@r2 start_time=1970-01-01T00:00:05Z du=3.0\nAAAAAA\n+\nIIIIII\n@r3 start_time=1970-01-01T00:00:10Z du=4.0\nGGGGGGGG\n+\n55555555\n" > tiny.fastq && \
		$(RUN_FASTCAT) tiny.fastq --output out --time-bin-size 10 --time-bin-k 64 > /dev/null && \
		test -f out/time_bin_summary.tsv && \
		test $$(wc -l < out/time_bin_summary.tsv) -eq 3 && \
		awk 'BEGIN { FS = OFS = "\t"; prev_reads = -1; prev_bases = -1; } \
			NR == 1 { if (NF != 24 || $$1 != "time_bin_start_s" || $$24 != "translocation_speed_q95") exit 1; next; } \
			{ \
				if (NF != 24) exit 1; \
				if ($$1 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$2 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/) exit 1; \
				if ($$3 < 0 || $$4 < 0 || $$5 < 0 || $$6 < 0) exit 1; \
				if ($$19 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$20 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$21 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$22 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$23 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $$24 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/) exit 1; \
				if (prev_reads >= 0 && ($$5 < prev_reads || $$6 < prev_bases)) exit 1; \
				prev_reads = $$5; prev_bases = $$6; \
				reads += $$3; bases += $$4; rows++; \
			} \
			END { \
				if (rows != 2) exit 1; \
				if (reads != 3) exit 1; \
				if (bases != 18) exit 1; \
				if (prev_reads != 3) exit 1; \
				if (prev_bases != 18) exit 1; \
			}' out/time_bin_summary.tsv && \
		printf "@u start_time=1970-01-01T00:00:00Z du=2.0\nACGT\n+\nIIII\n@b barcode=barcode0001 start_time=1970-01-01T00:00:05Z du=4.0\nAAAA\n+\nIIII\n" > demux.fastq && \
		$(RUN_FASTCAT) demux.fastq --output demux_out --demultiplex --time-bin-size 10 --time-bin-k 64 > /dev/null && \
		test -f demux_out/unclassified/time_bin_summary.tsv && \
		test -f demux_out/barcode0001/time_bin_summary.tsv)
test_fastcat_output: test_fastcat_time_bin_summary


.PHONY: test_fastcat_channel_summary
test_fastcat_channel_summary: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat channel_summary output shape and totals)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/channel_summary, \
		printf "@r1 ch=12 mux=1\nACGT\n+\nIIII\n@r2 ch=12 mux=1\nAAAAAA\n+\nIIIIII\n@r3 ch=12 mux=3\nGGG\n+\nIII\n@r4 ch=0 mux=2\nTT\n+\nII\n@r5 ch=7 mux=0\nCC\n+\nII\n" > tiny.fastq && \
		$(RUN_FASTCAT) tiny.fastq --output out > /dev/null && \
		test -f out/channel_summary.tsv && \
		awk 'BEGIN { FS = OFS = "\t"; rows = 0; reads = 0; bases = 0; } \
			NR == 1 { if (NF != 4 || $$1 != "channel" || $$4 != "base_count") exit 1; next; } \
			{ \
				if (NF != 4) exit 1; \
				if ($$1 < prev_channel || ($$1 == prev_channel && $$2 < prev_mux)) exit 1; \
				prev_channel = $$1; prev_mux = $$2; \
				rows++; reads += $$3; bases += $$4; \
				if ($$1 == 7 && $$2 == 0 && $$3 != 1) exit 1; \
				if ($$1 == 7 && $$2 == 0 && $$4 != 2) exit 1; \
				if ($$1 == 12 && $$2 == 1 && $$3 != 2) exit 1; \
				if ($$1 == 12 && $$2 == 1 && $$4 != 10) exit 1; \
				if ($$1 == 12 && $$2 == 3 && $$3 != 1) exit 1; \
				if ($$1 == 12 && $$2 == 3 && $$4 != 3) exit 1; \
			} \
			END { if (rows != 3 || reads != 4 || bases != 15) exit 1; }' out/channel_summary.tsv && \
		printf "@u ch=1 mux=1\nACGT\n+\nIIII\n@b barcode=barcode0001 ch=2 mux=4\nAAAA\n+\nIIII\n" > demux.fastq && \
		$(RUN_FASTCAT) demux.fastq --output demux_out --demultiplex > /dev/null && \
		test -f demux_out/unclassified/channel_summary.tsv && \
		test -f demux_out/barcode0001/channel_summary.tsv)
test_fastcat_output: test_fastcat_channel_summary


.PHONY: test_fastcat_summary_json
test_fastcat_summary_json: $(FASTCAT_BIN) $(REPORT_VENV)/.deps test_summary_json_validator
	@$(call TEST_DESC,Fastcat summary.json structure and totals)
	rm -rf $(TEST_TMP_DIR)/summary_json_fixture; \
	mkdir -p $(TEST_TMP_DIR)/summary_json_fixture; \
	pushd $(TEST_TMP_DIR)/summary_json_fixture > /dev/null; \
	printf "@r1 start_time=1970-01-01T00:00:00Z ch=12 mux=1 du=1.0\nACGT\n+\n!!!!\n@r2 barcode=barcode0001 start_time=1970-01-01T00:00:05Z ch=12 mux=3 du=2.0\nAAAAAA\n+\nIIIIII\n@r3 barcode=barcode0001 start_time=1970-01-01T00:00:15Z ch=9 mux=2 du=4.0\nGGGGGGGG\n+\n55555555\n" > tiny.fastq; \
	$(RUN_FASTCAT) $(TEST_TMP_DIR)/summary_json_fixture/tiny.fastq --output $(TEST_TMP_DIR)/summary_json_fixture/out --demultiplex --time-bin-size 10 --time-bin-k 64 --json-summary > /dev/null; \
	test -f $(TEST_TMP_DIR)/summary_json_fixture/out/unclassified/summary.json; \
	test -f $(TEST_TMP_DIR)/summary_json_fixture/out/barcode0001/summary.json; \
	$(call COMPARE_JSON,$(TEST_TMP_DIR)/summary_json_fixture/out/barcode0001/summary.json,$(TEST_EXPECTED_SUMMARY_JSON)/summary_json.barcode0001.json.gz) && \
	$(call COMPARE_JSON,$(TEST_TMP_DIR)/summary_json_fixture/out/unclassified/summary.json,$(TEST_EXPECTED_SUMMARY_JSON)/summary_json.unclassified.json.gz); \
	popd > /dev/null; \
	rm -rf $(TEST_TMP_DIR)/summary_json_fixture
test_fastcat_output: test_fastcat_summary_json


.PHONY: test_summary_json_validator
test_summary_json_validator: $(REPORT_VENV)/.deps tools/workflow_glue/validate_summary_json.py $(SUMMARY_JSON_SCHEMA)
	@$(call TEST_DESC,Summary JSON fixtures validate)
	$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) $(TEST_EXPECTED_SUMMARY_JSON)/summary_json.barcode0001.json.gz $(TEST_EXPECTED_SUMMARY_JSON)/summary_json.unclassified.json.gz


.PHONY: test_fastcat_summary_json_counters
test_fastcat_summary_json_counters: $(FASTCAT_BIN) $(REPORT_VENV)/.deps tools/workflow_glue/validate_summary_json.py $(SUMMARY_JSON_SCHEMA)
	@$(call TEST_DESC,Fastcat summary.json counter arrays)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/summary_json_counters, \
		printf "@r1 runid=12345678-1234-1234-1234-123456789012 basecall_model_version_id=model-a ch=1 mux=1 start_time=1970-01-01T00:00:00Z du=1.0\\nACGT\\n+\\nIIII\\n" > reads.fastq && \
		$(RUN_FASTCAT) reads.fastq --output out --json-summary > /dev/null && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) out/summary.json && \
		$(REPORT_PY) -c 'import json; data = json.load(open("out/summary.json")); assert len(data["counts"]["run_ids"]) == 1; assert len(data["counts"]["basecallers"]) == 1')
test_fastcat_output: test_fastcat_summary_json_counters


.PHONY: test_fastcat_regression_results
test_fastcat_regression_results: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat regression results)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/regression, \
		ln -s $(TEST_FIXTURES_FASTQ)/basic data && mkdir -p work && \
		pushd work > /dev/null && \
		$(RUN_FASTCAT) ../data/*.fastq.gz --sample sample \
			--output $(TEST_TMP_DIR)/fastcat/regression/out1 \
			--per-read \
			> $(TEST_TMP_DIR)/fastcat/regression/concat.fastq && \
		$(RUN_FASTCAT) ../data/*.fastq.gz \
			--output $(TEST_TMP_DIR)/fastcat/regression/out2 \
			> $(TEST_TMP_DIR)/fastcat/regression/concat.reheader.fastq && \
		popd > /dev/null && \
		$(call COMPARE_TABLES,$(TEST_TMP_DIR)/fastcat/regression/out1/per_read.tsv,$(TEST_EXPECTED_FASTCAT)/per-read-stats.tsv,per-read) && \
		$(call COMPARE_TABLES,$(TEST_TMP_DIR)/fastcat/regression/out1/per_file.tsv,$(TEST_EXPECTED_FASTCAT)/per-file-stats.tsv,per-file) && \
		$(call COMPARE_FASTQ,concat.fastq,$(TEST_EXPECTED_FASTCAT)/concat.sorted.fastq.gz,concat) && \
		$(call COMPARE_FASTQ,concat.reheader.fastq,$(TEST_EXPECTED_FASTCAT)/concat.reheader.sorted.fastq.gz,concat.reheader))
test_fastcat_equivalence: test_fastcat_regression_results


.PHONY: test_fastcat_perfile_all_filtered
test_fastcat_perfile_all_filtered: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat per-file stats when all reads filtered)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/perfile_filtered, \
		printf "@r1\nACGT\n+\n!!!!\n@r2\nACGT\n+\n!!!!\n" > short.fastq && \
		$(RUN_FASTCAT) short.fastq --min-length 100 --sample sample --output out > /dev/null && \
		awk 'NR==1{for(i=1;i<=NF;i++)h[$$i]=i;next} NR==2{ \
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_seen"]!=2||$$h["r_record_ok"]!=0||$$h["r_too_short"]!=2){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)
test_fastcat_filters: test_fastcat_perfile_all_filtered


.PHONY: test_fastcat_stdin_data_equiv
test_fastcat_stdin_data_equiv: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat stdin data matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/stdin_data_equiv, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output out_file > out_file.fastq && \
		cat $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz | $(RUN_FASTCAT) - --output out_stdin > out_stdin.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_stdin.fastq,fc_stdin_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,fc_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,fc_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,fc_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)
test_fastcat_inputs: test_fastcat_stdin_data_equiv


.PHONY: test_fastcat_fofn_equiv
test_fastcat_fofn_equiv: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat --fofn matches explicit file arguments)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/fofn_equiv, \
		printf "%s\n%s\n" $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz > paths.txt && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz --output out_file > out_file.fastq && \
		cat paths.txt | $(RUN_FASTCAT) --fofn - --output out_fofn > out_fofn.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_fofn.fastq,fc_fofn_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_fofn/runids.tsv,fc_fofn_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_fofn/basecallers.tsv,fc_fofn_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_fofn/per_file.tsv,fc_fofn_perfile) && \
		diff -r out_file/histograms out_fofn/histograms)
test_fastcat_inputs: test_fastcat_fofn_equiv


.PHONY: test_fastcat_stdin_requires_dash
test_fastcat_stdin_requires_dash: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat stdin modes require explicit '-')
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/stdin_requires_dash, \
		set +e; \
		cat $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz | $(RUN_FASTCAT) --output out_missing_dash > /dev/null 2>&1; \
		rc_data=$$?; \
		printf "%s\n%s\n" $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz | \
			$(RUN_FASTCAT) --fofn --output out_missing_fofn_dash > /dev/null 2>&1; \
		rc_fofn=$$?; \
		set -e; \
		test $$rc_data -ne 0 && test $$rc_fofn -ne 0)
test_fastcat_inputs: test_fastcat_stdin_requires_dash


.PHONY: test_fastcat_nested_output_dir
test_fastcat_nested_output_dir: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat supports nested output directories)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/nested_output, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output nested/out/path > /dev/null && \
		test -f nested/out/path/per_file.tsv && \
		test -f nested/out/path/histograms/length.hist)
test_fastcat_output: test_fastcat_nested_output_dir


.PHONY: test_fastcat_bam_input_equiv
test_fastcat_bam_input_equiv: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat BAM input (unaligned) matches FASTQ input)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_input_equiv, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz --output out_fastq > out_fastq.fastq && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz --output tmp2 -B > b2.bam && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam b2.bam --output out_bam > out_bam.fastq && \
		cat out_fastq.fastq \
			| sed -E 's/\t[qQ][sS]:f:[^\t[:space:]]+//g' \
			| $(DEV_TOOLS) sort_fastq > fc_bam_input_concat.sorted && \
		cat out_bam.fastq \
			| sed -E 's/\t[qQ][sS]:f:[^\t[:space:]]+//g' \
			| $(DEV_TOOLS) sort_fastq > fc_bam_input_concat.expected.sorted && \
		diff fc_bam_input_concat.sorted fc_bam_input_concat.expected.sorted && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/runids.tsv,out_bam/runids.tsv,fc_bam_input_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/basecallers.tsv,out_bam/basecallers.tsv,fc_bam_input_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_fastq/per_file.tsv,out_bam/per_file.tsv,fc_bam_input_perfile) && \
		diff -r out_fastq/histograms out_bam/histograms)
test_fastcat_inputs: test_fastcat_bam_input_equiv


.PHONY: test_fastcat_bam_input_fofn_equiv
test_fastcat_bam_input_fofn_equiv: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat BAM input via --fofn matches explicit BAM args)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_input_fofn_equiv, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc2.fastq.gz --output tmp2 -B > b2.bam && \
		printf "%s\n%s\n" b1.bam b2.bam > bam_paths.txt && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam b2.bam --output out_file > out_file.fastq && \
		$(RUN_FASTCAT) --input-fmt bam --fofn bam_paths.txt --output out_fofn > out_fofn.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_fofn.fastq,fc_bam_fofn_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_fofn/runids.tsv,fc_bam_fofn_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_fofn/basecallers.tsv,fc_bam_fofn_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_fofn/per_file.tsv,fc_bam_fofn_perfile) && \
		diff -r out_file/histograms out_fofn/histograms)
test_fastcat_inputs: test_fastcat_bam_input_fofn_equiv


.PHONY: test_fastcat_bam_input_stdin_equiv
test_fastcat_bam_input_stdin_equiv: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat BAM input from stdin matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_input_stdin_equiv, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output tmp1 -B > b1.bam && \
		$(RUN_FASTCAT) --input-fmt bam b1.bam --output out_file > out_file.fastq && \
		cat b1.bam | $(RUN_FASTCAT) --input-fmt bam - --output out_stdin > out_stdin.fastq && \
		$(call COMPARE_FASTQ,out_file.fastq,out_stdin.fastq,fc_bam_stdin_concat) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,fc_bam_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,fc_bam_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,fc_bam_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)
test_fastcat_inputs: test_fastcat_bam_input_stdin_equiv


.PHONY: test_fastcat_bam_input_rejects_mapped
test_fastcat_bam_input_rejects_mapped: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat BAM input rejects mapped BAM)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_input_rejects_mapped, \
		set +e; \
		$(RUN_FASTCAT) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --output out_bad > out.log 2>&1; \
		rc=$$?; \
		set -e; \
		test $$rc -ne 0 && grep -q "only accepts unaligned BAM input" out.log)
test_fastcat_inputs: test_fastcat_bam_input_rejects_mapped


.PHONY: mem_check_fastcat
mem_check_fastcat: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck fastcat (fastq))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/memcheck/fastq, \
		$(GRIND) $(FASTCAT_BIN) fastq $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz > /dev/null)
test_fastcat_memory: mem_check_fastcat


.PHONY: mem_check_fastcat_bam
mem_check_fastcat_bam: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck fastcat (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/memcheck/bam, \
		$(GRIND) $(FASTCAT_BIN) fastq $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz -B > /dev/null)
test_fastcat_memory: mem_check_fastcat_bam


.PHONY: mem_check_fastcat_demultiplex
mem_check_fastcat_demultiplex: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck fastcat demultiplex (fastq))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/memcheck/demux, \
		$(GRIND) $(FASTCAT_BIN) fastq $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz --demultiplex --output demultiplex > /dev/null)
test_fastcat_memory: mem_check_fastcat_demultiplex


.PHONY: mem_check_fastcat_demultiplex_bam
mem_check_fastcat_demultiplex_bam: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck fastcat demultiplex (bam))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/memcheck/demux_bam, \
		$(GRIND) $(FASTCAT_BIN) fastq $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz --demultiplex --output demultiplex -B > /dev/null)
test_fastcat_memory: mem_check_fastcat_demultiplex_bam


###
# bamstats tests

.PHONY: test_bamstats_demux_single_equiv
test_bamstats_demux_single_equiv: $(FASTCAT_BIN) $(REPORT_VENV)/.deps
	@$(call TEST_DESC,Bamstats demux vs non-demux (single barcode))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/demux_single, \
		$(call TEST_DEMUX_SINGLE_EQUIV,\
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output tmp_fc --per-read -B > bc1.bam && \
			$(RUN_BAMSTATS) bc1.bam --per-read --output plain_out --json-summary > /dev/null && \
			$(RUN_BAMSTATS) bc1.bam --per-read --demultiplex --output demux_out --json-summary > /dev/null,\
			$(BAMSTATS_EXPECTED_FILES),\
			bs_demux_single,\
			COMPARE_TABLES))
test_bamstats_core: test_bamstats_demux_single_equiv


.PHONY: test_meta_bamstats
test_meta_bamstats: $(FASTCAT_BIN)
	@$(call TEST_DESC,Metadata bamstats)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/meta/bam, \
		$(call TEST_METADATA,$(RUN_BAMSTATS),$(TEST_FIXTURES_BAM)/metadata,parse_rg,*.bam,$(TEST_EXPECTED_BAMSTATS)/metadata))
test_meta: test_meta_bamstats


.PHONY: test_bamstats_sample_column
test_bamstats_sample_column: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats sample_name column only when --sample set)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/sample_column, \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --output out_nosample > /dev/null && \
		! head -n 1 out_nosample/flagstats.tsv | grep -q "sample_name" && \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --sample "" --output out_empty > /dev/null && \
		head -n 1 out_empty/flagstats.tsv | grep -q "sample_name" && \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --sample foo --output out_sample > /dev/null && \
		head -n 1 out_sample/flagstats.tsv | grep -q "sample_name")
test_bamstats_core: test_bamstats_sample_column


.PHONY: test_bamstats_summary_json_schema
test_bamstats_summary_json_schema: $(FASTCAT_BIN) $(REPORT_VENV)/.deps tools/workflow_glue/validate_summary_json.py $(SUMMARY_JSON_SCHEMA)
	@$(call TEST_DESC,Bamstats summary.json schema including flagstats)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/summary_json_schema, \
		printf "@r1 runid=12345678-1234-1234-1234-123456789012 basecall_model_version_id=model-a ch=1 mux=1 start_time=1970-01-01T00:00:00Z du=1.0\\nACGT\\n+\\nIIII\\n" > reads.fastq && \
		$(RUN_FASTCAT) reads.fastq --output fastcat_out -B > reads.bam && \
		$(RUN_BAMSTATS) reads.bam --output out --json-summary > /dev/null && \
		test -s out/summary.json && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) out/summary.json && \
		$(REPORT_PY) -c 'import json; assert json.load(open("out/summary.json"))["flagstats"]')
test_bamstats_core: test_bamstats_summary_json_schema


.PHONY: test_bamstats_coverage_demux_layout
test_bamstats_coverage_demux_layout: $(FASTCAT_BIN) samtools
	@$(call TEST_DESC,Bamstats coverage demux layout)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/coverage_demux, \
		cp $(TEST_FIXTURES_BAM)/coverage/coverage_demux.sam cov.sam && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --coverage --demultiplex --output out > /dev/null && \
		if [ -d out/coverage ]; then echo "Unexpected top-level coverage directory"; exit 1; fi && \
		test -f out/barcode0001/coverage/global.bed.gz && \
		test -f out/barcode0002/coverage/global.bed.gz && \
		test ! -d out/unclassified/coverage)
test_bamstats_coverage: test_bamstats_coverage_demux_layout


.PHONY: test_bamstats_stdin_equiv
test_bamstats_stdin_equiv: $(FASTCAT_BIN) samtools
	@$(call TEST_DESC,Bamstats stdin BAM matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/stdin_equiv, \
		cp $(TEST_FIXTURES_BAM)/coverage/coverage_demux.sam cov.sam && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --output out_file > /dev/null && \
		cat cov.bam | $(RUN_BAMSTATS) - --output out_stdin > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,bs_stdin_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,bs_stdin_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,bs_stdin_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)
test_bamstats_inputs: test_bamstats_stdin_equiv


.PHONY: test_bamstats_stdin_bed_equiv
test_bamstats_stdin_bed_equiv: $(FASTCAT_BIN) samtools
	@$(call TEST_DESC,Bamstats stdin with BED filter matches file input)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/stdin_bed_equiv, \
		cp $(TEST_FIXTURES_BAM)/coverage/coverage_demux.sam cov.sam && \
		printf "chr1\t0\t75\n" > regions.bed && \
		$(REPO_ROOT)/samtools view -bS cov.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --bed regions.bed --output out_file > /dev/null && \
		cat cov.bam | $(RUN_BAMSTATS) - --bed regions.bed --output out_stdin > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/runids.tsv,out_stdin/runids.tsv,bs_stdin_bed_runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/basecallers.tsv,out_stdin/basecallers.tsv,bs_stdin_bed_basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,out_file/per_file.tsv,out_stdin/per_file.tsv,bs_stdin_bed_perfile) && \
		diff -r out_file/histograms out_stdin/histograms)
test_bamstats_inputs: test_bamstats_stdin_bed_equiv


.PHONY: test_bamstats_NM
test_bamstats_NM: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats NM validation)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/nm, \
		set +e; \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/validation/implausible_nm.sam --output out_bad_implausible > /dev/null 2> err_implausible.log; \
		rc_implausible=$$?; \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/validation/inconsistent_nm.sam --output out_bad_inconsistent > /dev/null 2> err_inconsistent.log; \
		rc_inconsistent=$$?; \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/validation/test.sam --output out_zero > /dev/null 2> err_zero.log; \
		rc_zero=$$?; \
		set -e; \
		test $$rc_implausible -ne 0 && grep -Eq "appears to contain implausible alignment information|zero aligned match/length" err_implausible.log && \
		test $$rc_inconsistent -ne 0 && grep -q "inconsistent NM tag" err_inconsistent.log && \
		test $$rc_zero -eq 0)
test_bamstats_core: test_bamstats_NM


.PHONY: test_bamstats_no_nan_metrics
test_bamstats_no_nan_metrics: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats rejects alignments that produce non-finite metrics)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/nan_guard, \
		set +e; $(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/validation/nan_bad.sam --output out_bad_nan > /dev/null 2> err_nan.log; rc=$$?; set -e; \
		test $$rc -ne 0 && grep -q "zero aligned match/length" err_nan.log)
test_bamstats_core: test_bamstats_no_nan_metrics


.PHONY: test_bamstats_polya
test_bamstats_polya: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats polya)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/polya, \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/special/RCS-100A.bam --poly-a --output out_polya > /dev/null && \
		diff out_polya/histograms/polya.hist $(TEST_EXPECTED_BAMSTATS)/polya/RCS-100A.bam.polya.hist)
test_bamstats_core: test_bamstats_polya


.PHONY: test_bamstats_filters
test_bamstats_filters: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats filters)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/filters, \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/filters/bamstats_filters.sam --tag-name HP --tag-value 1 --per-read --output out_tag > /dev/null && \
		grep -q "^read1" out_tag/per_read.tsv && \
		! grep -q "^read2" out_tag/per_read.tsv && \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/filters/bamstats_filters.sam --read-group rg2 --per-read --output out_rg > /dev/null && \
		grep -q "^read2" out_rg/per_read.tsv && \
		! grep -q "^read1" out_rg/per_read.tsv && \
		$(RUN_BAMSTATS) $(TEST_FIXTURES_BAM)/filters/bamstats_filters.sam --tag-name HP --tag-value 1 --read-group rg1 --per-read --output out_both > /dev/null && \
		grep -q "^read1" out_both/per_read.tsv && \
		! grep -q "^read2" out_both/per_read.tsv)
test_bamstats_core: test_bamstats_filters


.PHONY: test_bamstats_vs_fastcat_unaligned
test_bamstats_vs_fastcat_unaligned: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats vs fastcat (unaligned))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/vs_fastcat, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bcMixed.fastq.gz --output fc_out --per-read -B > reads.bam && \
		$(RUN_BAMSTATS) reads.bam --per-read --output bs_out > /dev/null && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/runids.tsv,bs_out/runids.tsv,runids) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/basecallers.tsv,bs_out/basecallers.tsv,basecallers) && \
		$(call COMPARE_TABLES_SKIP_FIRST,fc_out/per_file.tsv,bs_out/per_file.tsv,per_file) && \
		cut -f1$(comma)3-8 fc_out/per_read.tsv > pr.fc && \
		cut -f1$(comma)3-8 bs_out/per_read.tsv > pr.bs && \
		$(call COMPARE_TABLES,pr.fc,pr.bs,per_read_non_demux) && \
		diff fc_out/histograms/length.hist bs_out/histograms/length.unmap.hist && \
		diff fc_out/histograms/quality.hist bs_out/histograms/quality.unmap.hist)
test_bamstats_equivalence: test_bamstats_vs_fastcat_unaligned


.PHONY: test_bamstats_vs_fastcat_unaligned_demux
test_bamstats_vs_fastcat_unaligned_demux: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats vs fastcat (unaligned, demux))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/vs_fastcat_demux, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bcMixed.fastq.gz --output fc_out --per-read --demultiplex > /dev/null && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bcMixed.fastq.gz --per-read -B > reads.bam && \
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
test_bamstats_equivalence: test_bamstats_vs_fastcat_unaligned_demux


.PHONY: test_bamstats_perfile_all_filtered
test_bamstats_perfile_all_filtered: $(FASTCAT_BIN)
	@$(call TEST_DESC,Bamstats per-file stats when all reads filtered)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/perfile_filtered, \
		printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\n" > short.sam && \
		$(RUN_BAMSTATS) short.sam --tag-name ZZ --tag-value 1 --output out > /dev/null && \
		awk 'NR==1{for(i=1;i<=NF;i++)h[$$i]=i;next} NR==2{ \
			if($$h["n_seqs"]!=0||$$h["n_bases"]!=0||$$h["min_length"]!=0||$$h["max_length"]!=0||$$h["mean_quality"]!=0||$$h["f_file_ok"]!=1||$$h["f_unknown_error"]!=0||$$h["r_record_seen"]!=2||$$h["r_record_ok"]!=0||$$h["r_too_short"]!=0){ \
				print "unexpected per_file.tsv: "$$0; exit 1}}' out/per_file.tsv)
test_bamstats_core: test_bamstats_perfile_all_filtered


.PHONY: mem_check_bamstats_good
mem_check_bamstats_good: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamstats with good data)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/memcheck/good, \
		$(GRIND) $(FASTCAT_BIN) bam $(TEST_FIXTURES_BAM)/metadata/dna_r10.4.1_e8.2_400bps_hac@v4.3.0.bam > /dev/null)
test_bamstats_memory: mem_check_bamstats_good


.PHONY: mem_check_bamstats_bad
mem_check_bamstats_bad: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamstats with bad data)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/memcheck/bad, \
		$(GRIND) $(FASTCAT_BIN) bam $(TEST_FIXTURES_BAM)/metadata/bad-ones.bam > /dev/null)
test_bamstats_memory: mem_check_bamstats_bad


.PHONY: mem_check_bamstats_qcfail
mem_check_bamstats_qcfail: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamstats qcfail)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/memcheck/qcfail, \
		$(GRIND) $(FASTCAT_BIN) bam $(TEST_FIXTURES_BAM)/special/400ecoli-with-qcfail.bam > /dev/null)
test_bamstats_memory: mem_check_bamstats_qcfail


.PHONY: mem_check_bamstats_duplex
mem_check_bamstats_duplex: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamstats duplex)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/memcheck/duplex, \
		$(GRIND) $(FASTCAT_BIN) bam $(TEST_FIXTURES_BAM)/special/310dx.bam > /dev/null)
test_bamstats_memory: mem_check_bamstats_duplex


.PHONY: mem_check_bamstats_coverage
mem_check_bamstats_coverage: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamstats coverage)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/memcheck/coverage, \
		$(GRIND) $(FASTCAT_BIN) bam --coverage $(TEST_FIXTURES_BAM)/basic/400ecoli.bam > /dev/null)
test_bamstats_memory: mem_check_bamstats_coverage


###
# bamindex tests

.PHONY: mem_check_bamindex-build
mem_check_bamindex-build: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamindex build)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamindex/build, \
		cp $(TEST_FIXTURES_BAM)/basic/400.bam . && \
		$(GRIND) $(FASTCAT_BIN) index build 400.bam > /dev/null)
test_bamindex_memory: mem_check_bamindex-build


.PHONY: mem_check_bamindex-dump
mem_check_bamindex-dump: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamindex dump)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamindex/dump, \
		cp $(TEST_FIXTURES_BAM)/basic/400.bam . && \
		$(GRIND) $(FASTCAT_BIN) index build 400.bam > /dev/null && \
		$(GRIND) $(FASTCAT_BIN) index dump 400.bam.bci > /dev/null)
test_bamindex_memory: mem_check_bamindex-dump


.PHONY: mem_check_bamindex-fetch
mem_check_bamindex-fetch: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck bamindex fetch)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamindex/fetch, \
		cp $(TEST_FIXTURES_BAM)/basic/400.bam . && \
		$(GRIND) $(FASTCAT_BIN) index build 400.bam > /dev/null && \
		$(GRIND) $(FASTCAT_BIN) index fetch 400.bam --chunk 5 > /dev/null)
test_bamindex_memory: mem_check_bamindex-fetch


###
# fastlint tests

.PHONY: mem_check_fastlint
mem_check_fastlint: $(FASTCAT_BIN)
	@$(call TEST_DESC,Memcheck fastlint)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastlint/memcheck, \
		$(GRIND) $(FASTCAT_BIN) lint $(TEST_FIXTURES_FASTQ)/basic/*.fastq.gz > /dev/null)
test_fastlint_memory: mem_check_fastlint
