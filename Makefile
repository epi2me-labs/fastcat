OS := $(shell uname)
ifeq ($(OS), Darwin)
    # mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
    ARGP_INC ?= -I$(shell brew --prefix argp-standalone)/include
    CFLAGS ?= -fpic -O3 ${ARGP_INC}
else
    ARGP ?=
    ARGP_INC ?=
    CFLAGS ?= -fpic -msse3 -O3 ${ARGP_INC}
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


src/%.o: src/%.c
	$(CC) -Isrc -Ihtslib -c -pthread -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $^ -o $@


fastcat: src/fastcat/main.o src/fastcat/args.o src/fastcat/writer.o src/fastqcomments.o src/common.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lz -lm $(EXTRA_LIBS) \
		-o $@


bamstats: src/bamstats/main.o src/bamstats/args.o src/bamstats/readstats.o src/bamstats/bamiter.o src/fastqcomments.o src/common.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

bamindex: src/bamindex/main.o src/bamindex/build_main.o src/bamindex/fetch_main.o src/bamindex/dump_main.o src/bamindex/index.o $(STATIC_HTSLIB)
	$(CC) -Isrc -Ihtslib -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto $(EXTRA_LIBS) \
		-o $@

.PHONY: clean
clean:
	rm -rf fastcat bamstats bamindex src/fastcat/*.o src/bamstats/*.o src/bamindex/*.o src/*.o


.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean

.PHONY: mem_check_fastcat
mem_check_fastcat: fastcat
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./fastcat test/data/*.fastq.gz > /dev/null

.PHONY: mem_check_bamstats
mem_check_bamstats: bamstats
	$(VALGRIND) --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamstats test/bamstats/400ecoli.bam

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

.PHONY: regression_test_fastcat
regression_test_fastcat: fastcat
	if [ -d test/test-tmp ]; then rm -r test/test-tmp; fi
	mkdir test/test-tmp && \
	cd test/test-tmp && \
	../../fastcat ../data -s sample -H -f per-file-stats.tsv -r per-read-stats.tsv \
		| paste -d '|' - - - - | sort | tr '|' '\n' | gzip > concat.sorted.fastq.gz && \
	bash -c 'diff <(sort per-file-stats.tsv) \
		<(sort ../fastcat_expected_results/per-file-stats.tsv)' && \
	bash -c 'diff <(sort per-read-stats.tsv) \
		<(sort ../fastcat_expected_results/per-read-stats.tsv)' && \
	zdiff concat.sorted.fastq.gz ../fastcat_expected_results/concat.sorted.fastq.gz
	rm -r test/test-tmp
