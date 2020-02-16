#!/bin/bash
echo "Generating PDF..."
pdflatex -jobname=test_results test_results.tex
bibtex test_results
pdflatex -jobname=test_results test_results.tex
pdflatex -jobname=test_results test_results.tex
