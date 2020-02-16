#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=implementation_details implementation_details.tex
bibtex implementation_details
pdflatex -jobname=implementation_details implementation_details.tex
pdflatex -jobname=implementation_details implementation_details.tex
