#!/bin/bash

dirname=$1
let no_png=`ls $dirname/*.png | wc -l`
let no_png_m1=$no_png-1
convert -delay 20 -loop 0 `ls $dirname/*.png && ls -r $dirname/*.png | tail -n $no_png_m1` $dirname/dive.gif
