#!/bin/bash

NAME=fastcat

gcc -I src -o $NAME src/main.c src/args.c src/fastqcomments.c src/writer.c -lz -lm

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
