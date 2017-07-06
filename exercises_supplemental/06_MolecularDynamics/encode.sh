#!/bin/bash

#mencoder \
#        "mf://moviedata/*.png" -mf fps=15 \
#        -o ausgabevideo.avi \
#        -ovc lavc \
#        -lavcopts vcodec=mpeg4:vbitrate=8000

## create movie
#rm moviedata/*
mencoder mf://moviedata/*.png -mf w=800:h=600:fps=25:type=png -ovc copy -oac copy -lavcopts vcodec=mpeg4:vbitrate=8000 -o output.avi
