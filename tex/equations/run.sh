#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=equations equations.tex
bibtex equations
pdflatex -jobname=equations equations.tex
pdflatex -jobname=equations equations.tex
