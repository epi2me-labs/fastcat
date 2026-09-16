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
	@printf "  \033[36mdesktop\033[0m      Self-contained local FASTQ/BAM report bundle\n"
	@printf "\n"
	@printf "\033[1mBuild Options:\033[0m\n"
	@printf "  \033[36mUSE_DEFLATE=1\033[0m    Use libdeflate compression library\n"
	@printf "  \033[36mNOTHREADS=1\033[0m      Disable threading support\n"
	@printf "  \033[36mNOPURRFLATE=1\033[0m    Use zlib FASTQ decompression instead of Purrflate\n"
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
REPORT_TEMPLATE_VENDOR_STAMP := $(REPORT_DIR)/.template-vendor.stamp
REPORT_TEMPLATE_HEADER := src/core/fc_report_template_data.h
REPORT_TEMPLATE_SOURCE := src/core/fc_report_template_data.c
REPORT_EXAMPLE_N ?= 2000
REPORT_DEMO_N ?= 100000
REPORT_DEMO_UBAM_N ?= $(REPORT_DEMO_N)
REPORT_DEMO_DURATION_S ?= 172800
REPORT_DEMO_DIR = $(TEST_TMP_DIR)/report_html
REPORT_DEMO_FASTQ = $(REPORT_DEMO_DIR)/reads.fastq.gz
REPORT_DEMO_BAM = $(REPORT_DEMO_DIR)/aligned_reads.bam
DEV_TOOLS_SCRIPT := $(REPO_ROOT)/tools/dev-tools
DEV_TOOLS := python3 $(DEV_TOOLS_SCRIPT)
SIMULATE_READS := $(REPORT_PY) $(DEV_TOOLS_SCRIPT) simulate_reads
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
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@3)/lib -L$(shell brew --prefix curl)/lib
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

# Native MinGW/UCRT builds need argp-standalone and a reduced static HTSlib
# feature set. Keep this platform setup independent of Desktop packaging so a
# normal `make` builds and tests a Windows executable.
WINDOWS_BUILD := $(filter MINGW% MSYS% CYGWIN%,$(OS))

# initialize customizable variables if not set
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
HTS_CONF_ARGS ?=
HTSLIB_SYSTEM_LIBS ?= -llzma -lbz2 -lcurl -lcrypto
DESKTOP_TARGET ?= local
DESKTOP_VERSION := $(shell sed -n 's/^const char \*argp_program_version = "\([^"]*\)";.*/\1/p' src/version.c)
DESKTOP_STAGE := $(BUILD_DIR)/desktop-stage/fastcat-$(DESKTOP_VERSION)-$(DESKTOP_TARGET)
DESKTOP_ARCHIVE := $(BUILD_DIR)/desktop/fastcat-$(DESKTOP_VERSION)-$(DESKTOP_TARGET).tar.gz
DESKTOP_SHA256 := shasum -a 256
ifneq ($(findstring windows,$(DESKTOP_TARGET)),)
DESKTOP_EXECUTABLE := fastcat.exe
DESKTOP_SHA256 := sha256sum
else
DESKTOP_EXECUTABLE := fastcat
endif
PROFILER_LIBDIR ?=
ARGP ?=
ARGP_INC ?=

# Windows does not provide argp as part of the MinGW runtime. Build the
# vendored standalone implementation into the build directory for every
# native Windows build, including ordinary release and unit-test builds.
ifneq ($(WINDOWS_BUILD),)
    ARGP_STANDALONE_DIR := third_party/argp-standalone
    ARGP_STANDALONE_BUILD := $(BUILD_DIR)/argp-standalone
    ARGP := $(ARGP_STANDALONE_BUILD)/libargp.a
    ARGP_INC := -I$(ARGP_STANDALONE_DIR)
    ARGP_STANDALONE_PATCH := third_party/patches/argp-standalone-windows.patch
    ARGP_STANDALONE_CPPFLAGS := -DHAVE_DECL_FPUTS_UNLOCKED=0 \
        -DHAVE_DECL_FWRITE_UNLOCKED=0 -DHAVE_DECL_PUTC_UNLOCKED=0
    HTS_CONF_ARGS += --disable-libcurl --disable-plugins --disable-bz2 --disable-lzma --without-libdeflate
    HTSLIB_SYSTEM_LIBS := -lregex -ltre -lintl -liconv -lws2_32
    EXTRA_LDFLAGS += -static
endif

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
#   - desktop   : Optimized, self-contained local FASTQ/BAM report bundle
#
# USE_DEFLATE=1: Use libdeflate compression library with htslib
# NOTHREADS=1: Disable threading support
# NOPURRFLATE=1: Use the zlib FASTQ reader even if Purrflate is available

BUILD_MODE ?= release

# validate BUILD_MODE
VALID_BUILD_MODES := release debug sanitize profile valgrind desktop
ifeq ($(filter $(BUILD_MODE),$(VALID_BUILD_MODES)),)
    $(error Invalid BUILD_MODE: $(BUILD_MODE). Valid modes: $(VALID_BUILD_MODES))
endif

# disable threading
NOTHREADS ?=
ifeq ($(NOTHREADS), 1)
    CFLAGS += -DNOTHREADS
endif

# Purrflate is an optional private submodule. Source users who clone Fastcat
# from GitHub will normally not have it, so use the established zlib/kseq
# reader automatically in that case. Internal builds retain Purrflate by
# default, and NOPURRFLATE=1 remains useful for explicit fallback testing.
NOPURRFLATE ?=
PURRFLATE_ENABLED := 0
ifneq ($(NOPURRFLATE),1)
    ifneq ($(wildcard third_party/purrflate/src/zngp_reader.c),)
        PURRFLATE_ENABLED := 1
    endif
endif
ifeq ($(PURRFLATE_ENABLED),0)
    CFLAGS += -DNOPURRFLATE
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

else ifeq ($(BUILD_MODE),desktop)
    # Desktop build: optimized local report tool, with no optional HTSlib transports/codecs.
    CFLAGS += -O3 -DNDEBUG
ifeq ($(WINDOWS_BUILD),)
    HTS_CONF_ARGS += --disable-libcurl --disable-plugins --disable-bz2 --disable-lzma --without-libdeflate
    HTSLIB_SYSTEM_LIBS :=
endif
    ifneq ($(findstring linux,$(DESKTOP_TARGET)),)
        EXTRA_LDFLAGS += -static
    endif
    ifneq ($(findstring macos,$(DESKTOP_TARGET)),)
        EXTRA_LDFLAGS :=
    endif

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

ifeq ($(PURRFLATE_ENABLED),1)
ifneq ($(NOTHREADS),1)
    INCLUDES += -Ithird_party/purrflate/src -Ithird_party/purrflate/src/include
    PURRFLATE_SRCS := \
        third_party/purrflate/src/zngp_candidates.c \
        third_party/purrflate/src/zngp_chunk.c \
        third_party/purrflate/src/zngp_chunk_find.c \
        third_party/purrflate/src/zngp_deflate_header.c \
        third_party/purrflate/src/zngp_direct.c \
        third_party/purrflate/src/zngp_input.c \
        third_party/purrflate/src/zngp_decode.c \
        third_party/purrflate/src/zngp_pool.c \
        third_party/purrflate/src/zngp_reader.c \
        third_party/purrflate/src/zngp_stream_pages.c \
        third_party/purrflate/src/zngp_win.c \
        third_party/purrflate/src/zngp_symbolic.c
    PURRFLATE_OBJS := $(patsubst third_party/purrflate/src/%.c,$(OBJ_DIR)/purrflate/%.o,$(PURRFLATE_SRCS))
endif
endif

# core library object files
CORE_SRCS := \
	src/core/fc_common.c \
	src/core/fc_log.c \
	src/core/fc_stats.c \
	src/core/fc_quantile.c \
	src/core/fc_kll.c \
	src/core/fc_alignment_summary.c \
	src/core/fc_channel_summary.c \
	src/core/fc_channel_map.c \
	src/core/fc_channel_map_data.c \
	src/core/fc_metric_qscore_2d.c \
	src/core/fc_accuracy_qscore_2d.c \
	src/core/fc_length_qscore_2d.c \
	src/core/fc_run_summary.c \
	src/core/fc_speed_qscore_2d.c \
	src/core/fc_time.c \
	src/core/fc_time_bin.c \
	src/core/kh_int64_map.c \
	src/core/fc_stats_engine.c \
	src/core/fc_fastqcomments.c \
	src/core/fc_fastq_parser.c \
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
	src/cli/fastlint/args.c \
	src/cli/summary_json_schema_data.c

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

ifdef ARGP_STANDALONE_BUILD
ARGP_STANDALONE_SOURCES := $(shell find $(ARGP_STANDALONE_DIR) -type f)

$(ARGP_STANDALONE_BUILD)/libargp.a: $(ARGP_STANDALONE_SOURCES) $(ARGP_STANDALONE_PATCH)
	@echo Building bundled argp-standalone
	rm -rf $(ARGP_STANDALONE_BUILD)
	mkdir -p $(ARGP_STANDALONE_BUILD)
	cp -R $(ARGP_STANDALONE_DIR)/. $(ARGP_STANDALONE_BUILD)/
	cd $(ARGP_STANDALONE_BUILD) && patch -p1 < $(REPO_ROOT)/$(ARGP_STANDALONE_PATCH) \
		&& autoreconf --install \
		&& CPPFLAGS="$(ARGP_STANDALONE_CPPFLAGS)" ./configure \
		&& $(MAKE) libargp.a
endif

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

$(OBJ_DIR)/purrflate/%.o: third_party/purrflate/src/%.c third_party/zlib-ng/zlib.h $(BUILD_MAKEFILES) | $(OBJ_DIR) $(DEP_DIR)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@))
	$(CC) $(INCLUDES) -c -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@) $< -o $@

src/cli/summary_json_schema_data.c src/cli/summary_json_schema_data.h: test/fixtures/misc/summary_json.schema.json $(DEV_TOOLS_SCRIPT) tools/workflow_glue/generate_json_schema_c.py
	$(DEV_TOOLS) generate_json_schema_c --input $< --header src/cli/summary_json_schema_data.h --source src/cli/summary_json_schema_data.c

$(OBJ_DIR)/cli/main.o: src/cli/summary_json_schema_data.h

# version object
$(VERSION_OBJ): src/version.c third_party/zlib-ng/zlib.h $(BUILD_MAKEFILES) | $(OBJ_DIR) $(DEP_DIR)
	@mkdir -p $(dir $@)
	@mkdir -p $(dir $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@))
	$(CC) $(INCLUDES) -c -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@) $< -o $@

# dependencies on htslib patch for files that need it
$(OBJ_DIR)/core/fc_fastqcomments.o: | $(HTSLIB_PATCH_STAMP)
$(OBJ_DIR)/core/fc_read.o: | $(HTSLIB_PATCH_STAMP)

$(CORE_LIB): $(CORE_OBJS) $(PURRFLATE_OBJS) $(BUILD_MAKEFILES) | $(BUILD_DIR)
	rm -f $@
	$(AR) rcs $@ $(CORE_OBJS) $(PURRFLATE_OBJS)

# include dependency files
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(CORE_OBJS))
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(PURRFLATE_OBJS))
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(CLI_OBJS))
-include $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(VERSION_OBJ))

# oh, so thats what you wanted to do
$(FASTCAT_BIN): $(VERSION_OBJ) $(CLI_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) third_party/zlib-ng/libz.a $(ARGP) $(BUILD_MAKEFILES) | $(BUILD_DIR)
	@echo "\x1b[1;32mLinking $@\x1b[0m"
	$(CC) -pthread $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$(VERSION_OBJ) $(CLI_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) \
		third_party/zlib-ng/libz.a \
		$(ARGP) \
		-lm -lz -lpthread $(HTSLIB_SYSTEM_LIBS) $(EXTRA_LIBS) \
		-o $@

.PHONY: desktop_bundle
ifeq ($(BUILD_MODE),desktop)
desktop_bundle: $(FASTCAT_BIN)  ## assemble a Desktop bundle (set DESKTOP_TARGET)
	@rm -rf $(DESKTOP_STAGE)
	@mkdir -p $(DESKTOP_STAGE) $(dir $(DESKTOP_ARCHIVE))
	cp $(FASTCAT_BIN) $(DESKTOP_STAGE)/$(DESKTOP_EXECUTABLE)
	cp tools/desktop/nextflow_schema.json $(DESKTOP_STAGE)/nextflow_schema.json
	cp tools/desktop/README.md $(DESKTOP_STAGE)/README.md
	cp LICENSE $(DESKTOP_STAGE)/
	cd $(DESKTOP_STAGE) && $(DESKTOP_SHA256) $(DESKTOP_EXECUTABLE) nextflow_schema.json README.md LICENSE > SHA256SUMS
	tar -C $(dir $(DESKTOP_STAGE)) -czf $(DESKTOP_ARCHIVE) $(notdir $(DESKTOP_STAGE))
	@echo "Desktop bundle written to $(DESKTOP_ARCHIVE)"
else
desktop_bundle:
	@echo "Use BUILD_MODE=desktop DESKTOP_TARGET=<platform> make desktop_bundle" >&2
	@false
endif

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
	@find src -type f \( -name '*.c' -o -name '*.h' \) ! -path 'src/vendor/*' ! -path 'src/cli/summary_json_schema_data.c' ! -path 'src/core/fc_report_template_data.c' -print0 | xargs -0 $(CLANG_FORMAT) -i \
		-style='$(CLANG_FORMAT_STYLE)'

.PHONY: clang_format_staged
clang_format_staged:  ## run clang-format on staged source files
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@tmp=$$(mktemp); \
	trap 'rm -f "$$tmp"' EXIT; \
	git diff --cached --name-only --diff-filter=ACMR -z -- ':(glob)src/**/*.c' ':(glob)src/**/*.h' ':(exclude)src/vendor/**' ':(exclude)src/cli/summary_json_schema_data.c' ':(exclude)src/core/fc_report_template_data.c' > "$$tmp"; \
	if [ ! -s "$$tmp" ]; then \
		exit 0; \
	fi; \
	xargs -0 $(CLANG_FORMAT) -i -style='$(CLANG_FORMAT_STYLE)' < "$$tmp"; \
	xargs -0 git add -- < "$$tmp"

.PHONY: clang_format_check
clang_format_check:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || { echo "$(CLANG_FORMAT) not found"; exit 1; }
	@find src -type f \( -name '*.c' -o -name '*.h' \) ! -path 'src/vendor/*' ! -path 'src/cli/summary_json_schema_data.c' ! -path 'src/core/fc_report_template_data.c' -print0 | xargs -0 $(CLANG_FORMAT) --dry-run -Werror \
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

.PHONY: summary_json_schema_check
summary_json_schema_check:  ## verify the bundled summary JSON Schema asset is current
	@tmp=$$(mktemp -d); \
	trap 'rm -rf "$$tmp"' EXIT; \
	$(DEV_TOOLS) generate_json_schema_c \
		--input test/fixtures/misc/summary_json.schema.json \
		--header "$$tmp/summary_json_schema_data.h" \
		--source "$$tmp/summary_json_schema_data.c"; \
	diff -u src/cli/summary_json_schema_data.h "$$tmp/summary_json_schema_data.h"; \
	diff -u src/cli/summary_json_schema_data.c "$$tmp/summary_json_schema_data.c"

.PHONY: summary_json_schema_staged
summary_json_schema_staged:  ## regenerate and stage the schema asset when related files are staged
	@if ! git diff --cached --quiet -- \
		test/fixtures/misc/summary_json.schema.json \
		tools/workflow_glue/generate_json_schema_c.py \
		src/cli/summary_json_schema_data.h \
		src/cli/summary_json_schema_data.c; then \
		$(DEV_TOOLS) generate_json_schema_c \
			--input test/fixtures/misc/summary_json.schema.json \
			--header src/cli/summary_json_schema_data.h \
			--source src/cli/summary_json_schema_data.c; \
		git add -- src/cli/summary_json_schema_data.h src/cli/summary_json_schema_data.c; \
	fi

.PHONY: summary_json_docs
summary_json_docs:  ## regenerate the summary.json reference embedded in README.md
	$(DEV_TOOLS) generate_summary_json_docs \
		--schema test/fixtures/misc/summary_json.schema.json \
		--readme README.md

.PHONY: summary_json_docs_check
summary_json_docs_check:  ## verify the README summary.json reference is current
	$(DEV_TOOLS) generate_summary_json_docs \
		--schema test/fixtures/misc/summary_json.schema.json \
		--readme README.md \
		--check

.PHONY: summary_json_docs_staged
summary_json_docs_staged:  ## regenerate and stage README documentation when its source is staged
	@if ! git diff --cached --quiet -- \
		README.md \
		test/fixtures/misc/summary_json.schema.json \
		tools/workflow_glue/generate_summary_json_docs.py; then \
		$(MAKE) -s summary_json_docs; \
		git add -- README.md; \
	fi


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

$(REPORT_DIR)/template_example.fastq: $(REPORT_VENV)/.deps $(DEV_TOOLS_SCRIPT) tools/workflow_glue/simulate_reads.py | $(REPORT_DIR)
	$(SIMULATE_READS) \
		--n $(REPORT_EXAMPLE_N) \
		--shape 2.6 \
		--scale 30000 \
		--q1 10 \
		--q2 20 \
		--q1-weight 0.55 \
		--seed 42 \
		> $@

$(REPORT_DIR)/template_alignment.sam: $(REPORT_VENV)/.deps $(DEV_TOOLS_SCRIPT) tools/workflow_glue/simulate_reads.py | $(REPORT_DIR)
	$(SIMULATE_READS) \
		--output-format sam \
		--sam-coordinate-sort \
		--n 200 \
		--shape 2.6 \
		--scale 500 \
		--run-duration-s 172800 \
		--seed 43 \
		> $@

$(REPORT_DIR)/template_bam/summary.json: $(FASTCAT_BIN) $(REPORT_DIR)/template_alignment.sam test/fixtures/bam/coverage/report_simulated_regions.bed | $(REPORT_DIR)
	rm -rf $(REPORT_DIR)/template_bam
	$(FASTCAT_BIN) bam $(REPORT_DIR)/template_alignment.sam \
		--output $(REPORT_DIR)/template_bam \
		--coverage \
		--segments 100 \
		--coverage-beds test/fixtures/bam/coverage/report_simulated_regions.bed \
		--coverage-names regions \
		> /dev/null
	test -s $@

$(REPORT_DIR)/template_inputs/metadata.json: $(DEV_TOOLS_SCRIPT) tools/workflow_glue/write_wf_report_inputs.py | $(REPORT_DIR)
	$(DEV_TOOLS) write_wf_report_inputs --template --out-dir $(REPORT_DIR)/template_inputs

$(REPORT_DIR)/template_inputs/run_summary_overrides.json: $(REPORT_DIR)/template_inputs/metadata.json
	test -s $@

$(REPORT_DIR)/example_inputs/metadata.json: $(DEV_TOOLS_SCRIPT) tools/workflow_glue/write_wf_report_inputs.py | $(REPORT_DIR)
	$(DEV_TOOLS) write_wf_report_inputs --out-dir $(REPORT_DIR)/example_inputs

$(REPORT_DIR)/template_wf_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/template_bam/summary.json $(REPORT_DIR)/template_inputs/metadata.json $(REPORT_DIR)/template_inputs/run_summary_overrides.json $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/summary_json.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--summary-json $(REPORT_DIR)/template_bam/summary.json \
		--run-summary-overrides $(REPORT_DIR)/template_inputs/run_summary_overrides.json \
		--metadata $(REPORT_DIR)/template_inputs/metadata.json \
		--versions $(REPORT_DIR)/template_inputs/versions \
		--params $(REPORT_DIR)/template_inputs/params.json \
		--wf_version __FASTCAT_VERSION_HTML__

$(REPORT_DIR)/example_wf_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/template_bam/summary.json $(REPORT_DIR)/example_inputs/metadata.json $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/summary_json.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--summary-json $(REPORT_DIR)/template_bam/summary.json \
		--metadata $(REPORT_DIR)/example_inputs/metadata.json \
		--versions $(REPORT_DIR)/example_inputs/versions \
		--params $(REPORT_DIR)/example_inputs/params.json \
		--wf_version dev

$(REPORT_DIR)/time_plots_report.html: $(REPORT_VENV)/.deps $(REPORT_DIR)/template_bam/summary.json $(REPORT_DIR)/example_inputs/metadata.json $(DEV_TOOLS_SCRIPT) tools/workflow_glue/report.py tools/workflow_glue/summary_json.py tools/workflow_glue/bokeh_relational.py tools/workflow_glue/bokeh_matrix.py tools/workflow_glue/bokeh_boxplot.py tools/workflow_glue/bokeh_ribbon.py tools/workflow_glue/bokeh_seqsummary.py tools/workflow_glue/util.py | $(REPORT_MPLCONFIGDIR)
	MPLCONFIGDIR=$(REPORT_MPLCONFIGDIR) \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) report $@ \
		--summary-json $(REPORT_DIR)/template_bam/summary.json \
		--metadata $(REPORT_DIR)/example_inputs/metadata.json \
		--versions $(REPORT_DIR)/example_inputs/versions \
		--params $(REPORT_DIR)/example_inputs/params.json \
		--wf_version dev

.PHONY: report_time_plots_preview
report_time_plots_preview: $(REPORT_DIR)/time_plots_report.html  ## build Python-only report preview with time-bin plots

.PHONY: report_template
report_template: report_vendor_template  ## build the wf-template-derived HTML report template

.PHONY: report_vendor_template
report_vendor_template: $(REPORT_TEMPLATE_VENDOR_STAMP)  ## vendor report template into C assets

$(REPORT_TEMPLATE_VENDOR_STAMP): $(REPORT_DIR)/template_wf_report.html $(DEV_TOOLS_SCRIPT) tools/workflow_glue/vendor_report_template.py $(REPORT_TEMPLATE_HEADER) $(REPORT_TEMPLATE_SOURCE) | $(REPORT_DIR)
	$(DEV_TOOLS) vendor_report_template \
		--input $(REPORT_DIR)/template_wf_report.html \
		--template-output $(REPORT_DIR)/fastcat_report_template.html \
		--header $(REPORT_TEMPLATE_HEADER) \
		--source $(REPORT_TEMPLATE_SOURCE) \
		--symbol fc_report_template_html
	@touch $@

.PHONY: report_template_check
report_template_check: $(REPORT_DIR)/template_wf_report.html $(DEV_TOOLS_SCRIPT) tools/workflow_glue/vendor_report_template.py  ## verify the vendored report template assets are current
	@set -e; \
	tmp=$$(mktemp -d); \
	trap 'rm -rf "$$tmp"' EXIT; \
	$(DEV_TOOLS) vendor_report_template \
		--input $(REPORT_DIR)/template_wf_report.html \
		--header "$$tmp/fc_report_template_data.h" \
		--source "$$tmp/fc_report_template_data.c" \
		--verify-source $(REPORT_TEMPLATE_SOURCE) \
		--symbol fc_report_template_html; \
	diff -u $(REPORT_TEMPLATE_HEADER) "$$tmp/fc_report_template_data.h"

.PHONY: report_template_staged
report_template_staged:  ## refresh staged report-template C assets when their inputs change
	@if git diff --cached --quiet -- \
		Makefile \
		tools/requirements.txt \
		tools/workflow_glue/__init__.py \
		tools/workflow_glue/vendor_report_template.py \
		tools/workflow_glue/report.py \
		tools/workflow_glue/summary_json.py \
		tools/workflow_glue/bokeh_relational.py \
		tools/workflow_glue/bokeh_matrix.py \
		tools/workflow_glue/bokeh_boxplot.py \
		tools/workflow_glue/bokeh_ribbon.py \
		tools/workflow_glue/bokeh_seqsummary.py \
		tools/workflow_glue/simulate_reads.py \
		tools/workflow_glue/write_wf_report_inputs.py \
		tools/workflow_glue/util.py \
		test/fixtures/bam/coverage/report_simulated_regions.bed \
		src/core/fc_report_template_data.c \
		src/core/fc_report_template_data.h; then \
		exit 0; \
	fi; \
	$(MAKE) report_vendor_template; \
	git add -- src/core/fc_report_template_data.c src/core/fc_report_template_data.h

.PHONY: report_scaffold
report_scaffold: report_vendor_template  ## rebuild the developer HTML report scaffold


###
# Testing!
#

.PHONY: test
test: test_functional test_unit  ## run all tests

# unit tests use unity framework

UNIT_TEST_SRCS = test/unit/main.c $(sort $(wildcard test/unit/test_*.c)) test/unity/unity.c
ifeq ($(PURRFLATE_ENABLED),0)
UNIT_TEST_SRCS := $(filter-out test/unit/test_fastq_parser.c,$(UNIT_TEST_SRCS))
endif
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

test/unit/unit_tests: $(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) third_party/zlib-ng/libz.a $(ARGP) $(BUILD_MAKEFILES)
	$(CC) $(INCLUDES) $(WARNINGS) -fstack-protector-strong $(FORTIFY_CFLAG) \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$(UNIT_TEST_OBJS) $(CORE_LIB) $(STATIC_HTSLIB) third_party/zlib-ng/libz.a $(ARGP) \
		-lm -lz -lpthread $(HTSLIB_SYSTEM_LIBS) $(EXTRA_LIBS) \
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
.PHONY: test test_functional test_memory test_unit demo demo_prebuilt demo_ubam test_report test_report_missing_metadata
.PHONY: test_fastcat test_fastcat_memory test_fastcat_equivalence test_fastcat_demux test_fastcat_inputs test_fastcat_filters test_fastcat_output test_qc
.PHONY: test_bamstats test_bamstats_memory test_bamstats_core test_bamstats_equivalence test_bamstats_coverage test_bamstats_inputs
.PHONY: test_bamindex test_bamindex_memory
.PHONY: test_fastlint test_fastlint_memory
.PHONY: test_meta

.PHONY: test_functional
test_functional: test_fastcat test_bamstats test_bamindex test_fastlint test_meta test_qc  ## run integration tests

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
RUN_QC = $(PEPPER) $(FASTCAT_BIN) qc

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
FASTCAT_EXPECTED_FILES = runids.tsv basecallers.tsv per_read.tsv summary.json histograms/length.hist histograms/quality.hist
BAMSTATS_EXPECTED_FILES = runids.tsv basecallers.tsv per_read.tsv summary.json histograms/length.hist histograms/quality.hist histograms/accuracy.hist histograms/coverage.hist flagstats.tsv


###
# report tests

$(REPORT_DEMO_FASTQ): $(DEV_TOOLS_SCRIPT) tools/workflow_glue/simulate_reads.py $(REPORT_VENV)/.deps
	@mkdir -p $(@D)
	set -o pipefail; $(SIMULATE_READS) \
		--n $(REPORT_DEMO_N) \
		--run-duration-s $(REPORT_DEMO_DURATION_S) \
		--sample-id simulated-data \
		--max-len 1000000 \
		--seed 7 \
		| gzip -c > $@

$(REPORT_DEMO_BAM): $(DEV_TOOLS_SCRIPT) tools/workflow_glue/simulate_reads.py $(REPORT_VENV)/.deps
	@mkdir -p $(@D)
	$(SIMULATE_READS) \
		--output-format bam \
		--output-threads 4 \
		--sam-coordinate-sort \
		--n $(REPORT_DEMO_N) \
		--run-duration-s $(REPORT_DEMO_DURATION_S) \
		--sample-id simulated-alignment-data \
		--max-len 1000000 \
		--seed 17 \
		> $@

.PHONY: demo
demo: report_vendor_template  ## rebuild the current report template and generate HTML report demos
	$(MAKE) demo_prebuilt

.PHONY: demo_prebuilt
demo_prebuilt: $(FASTCAT_BIN) $(REPORT_DEMO_FASTQ) $(REPORT_DEMO_BAM) test/fixtures/bam/coverage/report_simulated_regions.bed  ## generate HTML report demos using the embedded template
	@$(call TEST_DESC,HTML report demo)
	rm -rf $(REPORT_DEMO_DIR)/fastcat_out $(REPORT_DEMO_DIR)/demux_out $(REPORT_DEMO_DIR)/aligned_out
	$(RUN_FASTCAT) $(REPORT_DEMO_FASTQ) \
		--output $(REPORT_DEMO_DIR)/fastcat_out \
		--html-report \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
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
	grep -q 'id="fastcat-alignment-section".*style="display:none;"' $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q 'id="fastcat-coverage-section".*style="display:none;"' $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "__FASTCAT_" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "fc_patch_fastcat_report" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Run configuration" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "Software versions" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "Workflow parameters" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "fastcat workflow provided by Oxford Nanopore Technologies" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	! grep -q "nextflow workflow" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	grep -q "Bokeh" $(TEST_TMP_DIR)/report_html/fastcat_out/fastcat_report.html
	$(RUN_FASTCAT) $(REPORT_DEMO_FASTQ) \
		--output $(REPORT_DEMO_DIR)/demux_out \
		--demultiplex \
		--html-report \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/fastcat_report.html
	test -s $(TEST_TMP_DIR)/report_html/demux_out/unclassified/summary.json
	$(RUN_BAMSTATS) $(REPORT_DEMO_BAM) \
		--output $(REPORT_DEMO_DIR)/aligned_out \
		--coverage \
		--coverage-beds test/fixtures/bam/coverage/report_simulated_regions.bed \
		--coverage-names regions \
		--html-report \
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q 'id="fastcat-alignment-section".*style=""' $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q 'id="fastcat-coverage-section".*style=""' $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -Eq 'data-fastcat-alignment-value="total-reads">[1-9][0-9]*<' $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -Eq 'data-fastcat-alignment-value="mapped-reads">[1-9][0-9]*<' $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -Eq 'data-fastcat-alignment-value="unmapped-reads">[1-9][0-9]*<' $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "Primary mappings" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "Alignment Concordance" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "Alignment concordance" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "Read alignment coverage" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "Whole reference" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	grep -q "User-provided: regions" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html
	! grep -q "__FASTCAT_" $(TEST_TMP_DIR)/report_html/aligned_out/fastcat_report.html

.PHONY: demo_ubam
demo_ubam: $(FASTCAT_BIN) $(REPORT_VENV)/.deps samtools  ## generate a synthetic unaligned BAM run and HTML report demo
	@$(call TEST_DESC,HTML report uBAM demo)
	rm -rf $(TEST_TMP_DIR)/report_html_ubam
	mkdir -p $(TEST_TMP_DIR)/report_html_ubam
	set -o pipefail; $(SIMULATE_READS) \
		--n $(REPORT_DEMO_UBAM_N) \
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
		> /dev/null
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	test -s $(TEST_TMP_DIR)/report_html_ubam/bam_out/summary.json
	$(REPORT_PY) -c 'import json; assert json.load(open("$(TEST_TMP_DIR)/report_html_ubam/bam_out/summary.json"))["run_summary"]["input_type"] == "unaligned bam"'
	grep -q "Input type" $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	grep -q "unaligned bam" $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	grep -q 'id="fastcat-alignment-section".*style="display:none;"' $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	grep -q 'id="fastcat-coverage-section".*style="display:none;"' $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	! grep -q "__FASTCAT_" $(TEST_TMP_DIR)/report_html_ubam/bam_out/fastcat_report.html
	rm -f $(TEST_TMP_DIR)/report_html_ubam/reads.fastq.gz

.PHONY: test_report
test_report: demo_prebuilt demo_ubam test_report_missing_metadata  ## run HTML report integration tests using the embedded template

.PHONY: test_report_missing_metadata
test_report_missing_metadata: $(FASTCAT_BIN) $(REPORT_VENV)/.deps
	@$(call TEST_DESC,HTML report null-state with missing metadata)
	rm -rf $(TEST_TMP_DIR)/report_missing_metadata
	mkdir -p $(TEST_TMP_DIR)/report_missing_metadata
	set -o pipefail; $(SIMULATE_READS) \
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
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output plain_out --per-read > /dev/null && \
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output demux_out --per-read --demultiplex > /dev/null,\
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
		test -f out_demux/unclassified/summary.json && \
		test -f out_demux/barcode0001/summary.json && \
		test -f out_demux/barcode0002/summary.json)
test_fastcat_demux: test_fastcat_demultiplex_outputs


.PHONY: test_fastcat_demultiplex_split
test_fastcat_demultiplex_split: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat demultiplex split files)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/demultiplex/split, \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --demultiplex --output out_split --reads-per-file 4 > /dev/null && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0000.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0000.fastq.gz,demux_split_0000) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0001.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0001.fastq.gz,demux_split_0001) && \
		$(call COMPARE_FASTQ,out_split/barcode0001/barcode0001_0002.fastq.gz,$(TEST_EXPECTED_DEMULTIPLEX)/barcode0001_0002.fastq.gz,demux_split_0002) && \
		test -f out_split/barcode0001/histograms/length.hist && \
		test -f out_split/barcode0001/summary.json)
test_fastcat_demux: test_fastcat_demultiplex_split


.PHONY: test_fastcat_summary_json
test_fastcat_summary_json: $(FASTCAT_BIN) $(REPORT_VENV)/.deps test_summary_json_validator
	@$(call TEST_DESC,Fastcat summary.json structure and totals)
	rm -rf $(TEST_TMP_DIR)/summary_json_fixture; \
	mkdir -p $(TEST_TMP_DIR)/summary_json_fixture; \
	pushd $(TEST_TMP_DIR)/summary_json_fixture > /dev/null; \
	printf "@r1 start_time=1970-01-01T00:00:00Z ch=12 mux=1 du=1.0\nACGT\n+\n!!!!\n@r2 barcode=barcode0001 start_time=1970-01-01T00:00:05Z ch=12 mux=3 du=2.0\nAAAAAA\n+\nIIIIII\n@r3 barcode=barcode0001 start_time=1970-01-01T00:00:15Z ch=9 mux=2 du=4.0\nGGGGGGGG\n+\n55555555\n" > tiny.fastq; \
	$(RUN_FASTCAT) $(TEST_TMP_DIR)/summary_json_fixture/tiny.fastq --output $(TEST_TMP_DIR)/summary_json_fixture/out --demultiplex --time-bin-size 10 --time-bin-k 64 > /dev/null; \
	test -f $(TEST_TMP_DIR)/summary_json_fixture/out/unclassified/summary.json; \
	test -f $(TEST_TMP_DIR)/summary_json_fixture/out/barcode0001/summary.json; \
	$(call COMPARE_JSON,$(TEST_TMP_DIR)/summary_json_fixture/out/barcode0001/summary.json,$(TEST_EXPECTED_SUMMARY_JSON)/summary_json.barcode0001.json.gz) && \
	$(call COMPARE_JSON,$(TEST_TMP_DIR)/summary_json_fixture/out/unclassified/summary.json,$(TEST_EXPECTED_SUMMARY_JSON)/summary_json.unclassified.json.gz); \
	popd > /dev/null; \
	rm -rf $(TEST_TMP_DIR)/summary_json_fixture
test_fastcat_output: test_fastcat_summary_json


.PHONY: test_summary_json_validator
test_summary_json_validator: $(REPORT_VENV)/.deps tools/workflow_glue/summary_json.py tools/workflow_glue/validate_summary_json.py $(SUMMARY_JSON_SCHEMA)
	@$(call TEST_DESC,Summary JSON fixtures validate)
	$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) $(TEST_EXPECTED_SUMMARY_JSON)/summary_json.barcode0001.json.gz $(TEST_EXPECTED_SUMMARY_JSON)/summary_json.unclassified.json.gz
	$(REPORT_PY) -c 'from tools.workflow_glue import summary_json; data = summary_json.load("$(TEST_EXPECTED_SUMMARY_JSON)/summary_json.barcode0001.json.gz"); channels = summary_json.decode_dataframe(data["channels"]["data"]); assert {"channel", "mean_quality", "mean_concordance_qscore"} <= set(channels); assert summary_json.decode_ndarray(data["heatmaps"]["read_length_quality"]["data"]).ndim == 2'

.PHONY: test_fastcat_summary_json_counters
test_fastcat_summary_json_counters: $(FASTCAT_BIN) $(REPORT_VENV)/.deps tools/workflow_glue/validate_summary_json.py $(SUMMARY_JSON_SCHEMA)
	@$(call TEST_DESC,Fastcat summary.json counter arrays)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/summary_json_counters, \
		printf "@r1 runid=12345678-1234-1234-1234-123456789012 basecall_model_version_id=model-a ch=1 mux=1 start_time=1970-01-01T00:00:00Z du=1.0\\nACGT\\n+\\nIIII\\n" > reads.fastq && \
		$(RUN_FASTCAT) reads.fastq --output out > /dev/null && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) out/summary.json && \
		$(REPORT_PY) -c 'import json; data = json.load(open("out/summary.json")); assert len(data["counts"]["run_ids"]) == 1; assert len(data["counts"]["basecallers"]) == 1' && \
		printf "@r2\\nACGT\\n+\\nIIII\\n" > missing_metadata.fastq && \
		$(RUN_FASTCAT) missing_metadata.fastq --output missing_metadata > /dev/null && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) missing_metadata/summary.json && \
		$(REPORT_PY) -c 'import json; data = json.load(open("missing_metadata/summary.json")); file = data["files"][0]; assert data["run_summary"]["run_ids"] is None; assert data["run_summary"]["basecaller_models"] is None; assert file["run_id_counts"][0]["value"] == "" and file["run_id_counts"][0]["count"] == 1; assert file["basecaller_counts"][0]["value"] == "" and file["basecaller_counts"][0]["count"] == 1; assert open("missing_metadata/runids.tsv").read().splitlines()[1].split(chr(9))[1] == ""; assert open("missing_metadata/basecallers.tsv").read().splitlines()[1].split(chr(9))[1] == ""')
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/summary_json_filtered, \
		printf "@short barcode=barcode0001\\nAAAA\\n+\\nIIII\\n@lowq barcode=barcode0001\\nAAAAAA\\n+\\n!!!!!!\\n@pass barcode=barcode0001\\nAAAAAA\\n+\\nIIIIII\\n" > filtered.fastq && \
		$(RUN_FASTCAT) filtered.fastq --output filtered --demultiplex --min-length 5 --min-qscore 20 --html-report > /dev/null && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) filtered/barcode0001/summary.json && \
		$(REPORT_PY) -c 'import json; row = json.load(open("filtered/barcode0001/summary.json"))["run_summary"]; assert (row["records_seen"], row["records_passed_filters"], row["records_too_long"], row["records_too_short"], row["records_low_quality"], row["records_dust_masked"]) == (3, 1, 0, 1, 1, 0)' && \
		grep -A 1 'Records seen' filtered/barcode0001/fastcat_report.html | grep -q '>3<' && \
		grep -A 1 'Too short' filtered/barcode0001/fastcat_report.html | grep -q '>1<')
test_fastcat_output: test_fastcat_summary_json_counters


.PHONY: test_fastcat_summary_output_layout
test_fastcat_summary_output_layout: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat summary JSON output layout)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/summary_output_layout, \
		printf "@r1 barcode=barcode0001\\nACGT\\n+\\nIIII\\n" > reads.fastq && \
		$(RUN_FASTCAT) reads.fastq --output normal > /dev/null && \
		$(RUN_FASTCAT) reads.fastq --output demux --demultiplex > /dev/null && \
		$(RUN_FASTCAT) reads.fastq --output no_legacy --no-legacy-outputs > /dev/null && \
		for directory in normal demux/barcode0001 no_legacy; do \
			test -s $$directory/summary.json || exit 1; \
			for removed in run_summary.tsv channel_summary.tsv length_qscore_2d.tsv speed_qscore_2d.tsv accuracy_qscore_2d.tsv time_bin_summary.tsv; do \
				test ! -e $$directory/$$removed || exit 1; \
			done; \
		done)
test_fastcat_output: test_fastcat_summary_output_layout


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


.PHONY: test_fastcat_bam_input_redirects_qc
test_fastcat_bam_input_redirects_qc: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat BAM input is redirected to QC)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/fastcat/bam_input_redirects_qc, \
		set +e; \
		$(RUN_FASTCAT) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam > detected.log 2>&1; \
		detected_rc=$$?; \
		$(RUN_FASTCAT) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam > option.log 2>&1; \
		option_rc=$$?; \
		set -e; \
		test $$detected_rc -ne 0 && grep -q "use 'fastcat qc'" detected.log && \
		test $$option_rc -ne 0 && grep -q "fastcat qc --input-fmt bam" option.log)
test_fastcat_inputs: test_fastcat_bam_input_redirects_qc


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
# fastcat qc tests
.PHONY: test_qc_fastq_default_output
test_qc_fastq_default_output: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat qc FASTQ streams reads unless an output path is requested)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/qc/fastq_default, \
		$(RUN_QC) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output out > reads.fastq && \
		test -s reads.fastq && \
		test -s out/summary.json && \
		test ! -e out/per_file.tsv && \
		test ! -d out/histograms && \
		$(RUN_QC) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --reads-output saved --output saved > saved.stdout && \
		test ! -s saved.stdout && \
		test -s saved/saved.fastq.gz && \
		gzip -dc saved/saved.fastq.gz | grep -q '^@' && \
		$(RUN_QC) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --bam-out --reads-output saved --output saved_bam > saved_bam.stdout && \
		test ! -s saved_bam.stdout && test -s saved_bam/saved.bam && \
		$(RUN_QC) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --demultiplex --reads-output ignored.fastq.gz --output demultiplex > demultiplex.stdout && \
		test ! -s demultiplex.stdout && \
		test ! -e demultiplex/ignored.fastq.gz && \
		test -n "$$({ find demultiplex -name '*.fastq.gz' -type f; })" && \
		cp $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz reads.bam && \
		$(RUN_QC) reads.bam --output detected_from_content > detected.fastq && \
		test -s detected.fastq && test -s detected_from_content/summary.json && \
		$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output legacy > legacy.fastq && \
		test -s legacy.fastq && \
		$(RUN_QC) $(TEST_FIXTURES_FASTQ)/basic --output directory > directory.fastq && \
		test -s directory.fastq && test -s directory/summary.json)
test_qc: test_qc_fastq_default_output


.PHONY: test_qc_fastq_input_matrix
test_qc_fastq_input_matrix: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat qc supports plain and gzip FASTQ from files and stdin)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/qc/fastq_input_matrix, \
		gzip -dc $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz > reads.fastq && \
		expected_reads=$$(grep -c '^@' reads.fastq) && \
		failed=0 && \
		check_case() { \
			case_name=$$1; case_status=$$2; \
			actual_reads=$$(grep -c '^@' $$case_name.fastq); \
			if test $$case_status -ne 0 || test "$$actual_reads" != "$$expected_reads"; then \
				printf '%s: exit %s expected %s reads got %s\n' "$$case_name" "$$case_status" "$$expected_reads" "$$actual_reads" >&2; \
				failed=1; \
			fi; \
		}; \
		for threads in 1 2; do \
			$(RUN_QC) --input-fmt fastq --threads $$threads --output plain_file_$$threads reads.fastq > plain_file_$$threads.fastq; check_case plain_file_$$threads $$?; \
			set -o pipefail; cat reads.fastq | $(RUN_QC) --input-fmt fastq --threads $$threads --output plain_stdin_$$threads - > plain_stdin_$$threads.fastq; check_case plain_stdin_$$threads $$?; \
			$(RUN_QC) --input-fmt fastq --threads $$threads --output gzip_file_$$threads $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz > gzip_file_$$threads.fastq; check_case gzip_file_$$threads $$?; \
			set -o pipefail; cat $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz | $(RUN_QC) --input-fmt fastq --threads $$threads --output gzip_stdin_$$threads - > gzip_stdin_$$threads.fastq; check_case gzip_stdin_$$threads $$?; \
		done; \
		test $$failed -eq 0)
test_qc: test_qc_fastq_input_matrix


.PHONY: test_qc_reference_bam_policy
test_qc_reference_bam_policy: $(FASTCAT_BIN) $(REPORT_VENV)/.deps
	@$(call TEST_DESC,Fastcat qc reference BAM is analysis-only)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/qc/reference_bam, \
		$(RUN_QC) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --min-length 999999 --output out > reads.fastq 2> log && \
		test ! -s reads.fastq && \
		test -s out/summary.json && \
		$(REPORT_PY) -c 'import json; assert json.load(open("out/summary.json"))["coverage"]["targets"]' && \
		test ! -d out/coverage && \
		test ! -e out/flagstats.tsv && \
		grep -q "filtering is disabled" log && \
		$(RUN_QC) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --legacy-outputs --output legacy > /dev/null && \
		test -f legacy/flagstats.tsv && \
		test -f legacy/coverage/global.summary.txt && \
		test -f legacy/coverage/global.dist.txt && \
		test ! -e legacy/coverage/global.bed.gz && \
		$(RUN_QC) $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --coverage --output traces > /dev/null && \
		test -f traces/coverage/global.bed.gz && \
		test ! -e traces/coverage/global.summary.txt && \
		set +e; \
		$(RUN_QC) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --reads-output saved --output rejected > /dev/null 2> rejected.log; \
		rejected_rc=$$?; \
		set -e; \
		test $$rejected_rc -ne 0 && \
		grep -q -- '--reads-output is not supported' rejected.log && \
		test ! -e rejected && test ! -e rejected.tmp)
test_qc: test_qc_reference_bam_policy


.PHONY: test_qc_multi_bam_restrictions
test_qc_multi_bam_restrictions: $(FASTCAT_BIN)
	@$(call TEST_DESC,Fastcat qc validates multi-BAM restrictions)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/qc/multi_bam, \
		set +e; \
		$(RUN_QC) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam $(TEST_FIXTURES_BAM)/special/310dx.bam --output mismatched > /dev/null 2> mismatch.log; \
		mismatch_rc=$$?; \
		$(RUN_QC) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam --coverage --output coverage > /dev/null 2> coverage.log; \
		coverage_rc=$$?; \
		$(RUN_QC) --input-fmt bam $(TEST_FIXTURES_BAM)/basic/400ecoli.bam $(TEST_FIXTURES_BAM)/special/310dx.bam --output retained --keep-failed-output > /dev/null 2> retained.log; \
		retained_rc=$$?; \
		set -e; \
		test $$mismatch_rc -ne 0 && grep -q "Reference dictionary" mismatch.log && \
		test $$coverage_rc -ne 0 && grep -q "Coverage is supported only with one BAM input" coverage.log && \
		test $$retained_rc -ne 0 && grep -q "Keeping partial output" retained.log && \
		test ! -e mismatched && test ! -e mismatched.tmp && \
		test ! -e coverage && test ! -e coverage.tmp && \
		test ! -e retained && test -d retained.tmp)
test_qc: test_qc_multi_bam_restrictions


###
# bamstats tests

.PHONY: test_bamstats_demux_single_equiv
test_bamstats_demux_single_equiv: $(FASTCAT_BIN) $(REPORT_VENV)/.deps
	@$(call TEST_DESC,Bamstats demux vs non-demux (single barcode))
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/demux_single, \
		$(call TEST_DEMUX_SINGLE_EQUIV,\
			$(RUN_FASTCAT) $(TEST_FIXTURES_FASTQ)/basic/bc1.fastq.gz --output tmp_fc --per-read -B > bc1.bam && \
			$(RUN_BAMSTATS) bc1.bam --per-read --output plain_out > /dev/null && \
			$(RUN_BAMSTATS) bc1.bam --per-read --demultiplex --output demux_out > /dev/null,\
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
	@$(call TEST_DESC,Bamstats summary.json schema including alignment)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/summary_json_schema, \
		printf "@r1 runid=12345678-1234-1234-1234-123456789012 basecall_model_version_id=model-a ch=1 mux=1 start_time=1970-01-01T00:00:00Z du=1.0\\nACGT\\n+\\nIIII\\n" > reads.fastq && \
		$(RUN_FASTCAT) reads.fastq --output fastcat_out -B > reads.bam && \
		$(RUN_BAMSTATS) reads.bam --output out > /dev/null && \
		test -s out/summary.json && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) out/summary.json && \
		$(REPORT_PY) -c 'import json; assert json.load(open("out/summary.json"))["alignment"]["per_reference_counts"]')
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
		test -f out/barcode0001/coverage/global.summary.txt && \
		test -f out/barcode0002/coverage/global.bed.gz && \
		test ! -d out/unclassified/coverage)
test_bamstats_coverage: test_bamstats_coverage_demux_layout


.PHONY: test_bamstats_coverage_summary_json
test_bamstats_coverage_summary_json: $(FASTCAT_BIN) samtools $(REPORT_VENV)/.deps test_summary_json_validator
	@$(call TEST_DESC,Bamstats coverage summary.json tables)
	$(call RUN_IN_TMP,$(TEST_TMP_DIR)/bamstats/coverage_summary_json, \
		cp $(TEST_FIXTURES_BAM)/coverage/coverage_demux.sam cov.sam && \
		sed -e '/^@/!s/$$/\tCH:i:1\tMX:i:1\tDU:f:1.0\tST:Z:1970-01-01T00:00:00Z/' cov.sam > tagged.sam && \
		printf 'chr1\t0\t200\n' > regions.bed && \
		$(REPO_ROOT)/samtools view -bS tagged.sam | $(REPO_ROOT)/samtools sort -o cov.bam && \
		$(REPO_ROOT)/samtools index cov.bam && \
		$(RUN_BAMSTATS) cov.bam --segments 100 --coverage-beds regions.bed --coverage-names regions --output out > /dev/null && \
		test -s out/summary.json && \
		test ! -e out/coverage/global.summary.txt && \
		test ! -e out/coverage/global.bed.gz && \
		$(REPORT_PY) $(DEV_TOOLS_SCRIPT) validate_summary_json --schema $(SUMMARY_JSON_SCHEMA) out/summary.json && \
		grep -q '"name":"global"' out/summary.json && \
		grep -q '"name":"segments_100"' out/summary.json && \
		grep -q '"name":"regions"' out/summary.json)
test_bamstats_coverage: test_bamstats_coverage_summary_json


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
		diff fc_out/histograms/length.hist bs_out/histograms/length.hist && \
		diff fc_out/histograms/quality.hist bs_out/histograms/quality.hist)
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
			g=bs_out/$$rel; \
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
