#!/usr/bin/env bash

## USAGE: from folder containing the pngs, run this. Output is output.mp4

ffmpeg -framerate 20 -i dump%04d.dat.png -c:v libx264 -s:v 1228x180 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
