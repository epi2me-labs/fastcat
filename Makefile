OS := $(shell uname)
ifeq ($(OS), Darwin)
	# mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L$(shell brew --prefix openssl@1.1)/lib
    ARGP ?= $(shell brew --prefix argp-standalone)/lib/libargp.a
	ARGP_INC ?= -I$(shell brew --prefix argp-standalone)/include
else
    ARGP ?=
	ARGP_INC ?=
endif


CC ?= gcc
CFLAGS ?= -fpic -msse3 -O3 ${ARGP_INC}
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
default: fastcat bamstats


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


fastcat: src/fastcat/main.o src/fastcat/args.o src/fastcat/writer.o src/fastqcomments.o
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


.PHONY: clean
clean:
	rm -rf fastcat bamstats src/fastcat/*.o src/bamstats/*.o src/*.o


.PHONY: clean_htslib
clean_htslib:
	cd htslib && make clean

.PHONY: mem_check_fastcat
mem_check_fastcat: fastcat
	valgrind --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./fastcat test/data/*.fastq.gz > /dev/null

.PHONY: mem_check_bamstats
mem_check_bamstats: bamstats
	valgrind --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s \
		./bamstats test/bamstats/400ecoli.bam
