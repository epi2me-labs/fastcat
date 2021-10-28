OS := $(shell uname)
ifeq ($(OS), Darwin)
	# mainly for dev builds using homebrew things
    EXTRA_LDFLAGS ?= -L/usr/local/Cellar/openssl@1.1/1.1.1k/lib
    ARGP ?= /usr/local/Cellar/argp-standalone/1.3/lib/libargp.a
else
    ARGP ?=
endif


CC ?= gcc
CFLAGS ?= -fpic -msse3 -O3 
EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=

fastcat: src/main.c src/args.c src/fastqcomments.c src/writer.c
	$(CC) -Isrc -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 \
		$(CFLAGS) $(EXTRA_CFLAGS) $(EXTRA_LDFLAGS) \
		$^ $(ARGP) \
		-lz -lm $(EXTRA_LIBS) \
		-o $@

.PHONY: clean
clean:
	rm -rf fastcat

.PHONY: mem_check
mem_check: fastcat
	valgrind --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s ./fastcat test_data/data/*.fastq.gz > /dev/null
