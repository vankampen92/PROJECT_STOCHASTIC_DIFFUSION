#!/bin/bash

# Usage:
# Example 1:  
# :~$ bash build-rename-movie.sh [movie name.avi]

# Correcting first time frame:
mv pgplot.png pgplot.png_1

# Renaming the png files in the right ordr
bash rename-png-files-in-order.sh 0

# Making the *.avi movie file
ffmpeg -framerate 25 -i %03d.png output.avi

# Renaming the output file with the input name
mv output.avi $1

# Deleting frames
rm *.png
rm pgplot*
