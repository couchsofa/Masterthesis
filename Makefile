all: clean
	rm -f *.pdf
	pdflatex  main.tex
	biber main.bcf
	pdflatex  main.tex
	biber main.bcf
	pdflatex  main.tex
	make clean
	chromium ./main.pdf

clean:
	rm -f *.aux *.log *.toc *.bbl *.run.xml *.bcf *.blg *.tdo *.out *.listing *.lof
