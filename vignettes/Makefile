all: MTBLS1582.html

%.html: %.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render("'$<'", "html_document")'
