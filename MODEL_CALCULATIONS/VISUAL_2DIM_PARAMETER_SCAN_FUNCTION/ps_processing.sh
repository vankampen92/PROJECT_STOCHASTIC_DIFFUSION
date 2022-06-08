rm output.*
ps2ps $1.ps output.ps
ps2eps output.ps
epstopdf output.eps
pdftoppm -png -rx 300 -ry 300 output.pdf output
mv output-1.png $2.png
