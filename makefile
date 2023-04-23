latex: datastep latex/main.tex
	cd latex && make

datastep: rmarkdown/document.Rmd
	rm /rmarkdown/document_cache/ -rf
	rm /rmarkdown/document_files/ -rf
	R CMD BATCH rmarkdown/knit.R
