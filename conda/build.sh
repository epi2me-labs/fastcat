#!/bin/bash

NAME=fastcat

gcc -o $NAME main.c -lz -lm

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
