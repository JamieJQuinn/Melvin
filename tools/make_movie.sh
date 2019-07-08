#!/usr/bin/env bash

## USAGE: from folder containing the pngs, run this. Output is output.mp4

cp ../../tools/movie_makefile ./makefile
make -j $(nproc)
