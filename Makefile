all: clean figures
	rm -f *.pdf
	pdflatex  main.tex
	biber main.bcf
	pdflatex  main.tex
	biber main.bcf
	pdflatex  main.tex
	make clean
	chromium ./main.pdf

clean:
	rm -f *.aux *.log *.toc *.bbl *.run.xml *.bcf *.blg *.tdo *.out *.listing *.lof *.tex.bak

figures:
	python plots/DampenerSweep.py
	python plots/Isolationsspektrum.py
	python plots/Isolationsspektrum2.py
	python plots/Isolationsspektrum3.py
	python plots/Antwortspektrum.py
	inkscape --export-dpi=300 -z -e images/voigt-kelvin-model.png svg/voigt-kelvin-model.svg
	inkscape --export-dpi=300 -z -e images/composition.png svg/composition.svg
	inkscape --export-dpi=300 -z -e images/composition2.png svg/composition2.svg
	inkscape --export-dpi=300 -z -e images/2MS_beispiel.png svg/2MS_beispiel.svg
	inkscape --export-dpi=300 -z -e images/1MS_AWS.png svg/1MS_AWS.svg