.PHONY: all, clean,

.SECONDARY:

all: paper.pdf

# This goes against all the nice rules below, but I couldn't get them to
# work.
paper.pdf: paper.tex paper.bib review-responses.tex
	pdflatex paper
	bibtex paper
	pdflatex paper

review-diff.tex: paper.tex
	latexdiff submitted_version.tex paper.tex > review-diff.tex

review-diff.pdf: review-diff.tex
	pdflatex review-diff.tex
	pdflatex review-diff.tex
	bibtex review-diff
	pdflatex review-diff.tex

clean:
	-rm -f *.aux *.log *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex -shell-escape $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex -interaction nonstopmode -shell-escape $<

%.bbl : %.aux paper.bib
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

