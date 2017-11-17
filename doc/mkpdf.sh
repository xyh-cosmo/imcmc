#!/bin/bash
TEX=xelatex
BIBFILE='EoS'
${TEX}  -synctex=1 ${BIBFILE}.tex
bibtex ${BIBFILE}.aux
${TEX}  -synctex=1 ${BIBFILE}.tex
${TEX}  -synctex=1 ${BIBFILE}.tex

