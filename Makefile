.PHONY: all, clean,

.SECONDARY:

all: paper.pdf



paper.pdf : paper.tex paper.bib

clean: 
	-rm *.aux *.log *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex -shell-escape $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex -interaction nonstopmode -shell-escape $<

%.bbl : %.aux bibs
	bibtex $<

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.svg
	inkscape $< --export-area-drawing --export-filename=$@
	# chromium --headless --no-pdf-header-footer --print-to-pdf=$@ $<
	# ./svg2pdf.sh $< $@

%.pdf : %.eps
	# inkscape $< --export-filename=$@
	epspdf $<

%.pdf : %.ink.svg
	inkscape $< --export-filename=$@

