#!/bin/bash

NAME=fastcat

gcc -o $NAME main.c args.c fastqcomments.c writer.c -lz -lm

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
