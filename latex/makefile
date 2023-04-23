run: main.tex
	pdflatex -halt-on-error --shell-escape main
	@bibtex main
	@pdflatex -halt-on-error --shell-escape main 2>&1 >/dev/null
	pdflatex -halt-on-error --shell-escape main
