OS := $(shell uname)

ifeq ($(OS), Darwin)
    ARGP ?= /usr/local/Cellar/argp-standalone/1.3/lib/libargp.a
else
    ARGP ?=
endif

fastcat: main.c args.c fastqcomments.c writer.c
	gcc -o fastcat $^ $(ARGP) -lz -lm


.PHONY: mem_check
mem_check: fastcat
	valgrind --error-exitcode=1 --tool=memcheck --leak-check=full --show-leak-kinds=all -s ./fastcat bc1.fastq.gz bc2.fastq.gz > /dev/null
