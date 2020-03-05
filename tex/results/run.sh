#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=hydro-mesh-results hydro-mesh-results.tex
bibtex hydro-mesh-results
pdflatex -jobname=hydro-mesh-results hydro-mesh-results.tex
pdflatex -jobname=hydro-mesh-results hydro-mesh-results.tex
