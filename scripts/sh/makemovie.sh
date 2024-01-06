#!/bin/bash

# make a movie out of png files in this directory using ffmpeg.

png_prefix=advection-1D-pwconst
video_fname=advection-piecewise-constant-1D-fixed-velocity.mp4

ffmpeg -framerate 12 -i $png_prefix-%04d-density-only.png -c:v libx264 -r 30 -pix_fmt yuv420p $video_fname
