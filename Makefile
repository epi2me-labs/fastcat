OS := $(shell uname)

ifeq ($(OS), Darwin)
    ARGP ?= /usr/local/Cellar//argp-standalone/1.3/lib/libargp.a
else
    ARGP ?=
endif

fastcat:
	gcc -o fastcat main.c $(ARGP) -lz -lm
