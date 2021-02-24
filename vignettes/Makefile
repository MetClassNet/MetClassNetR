all: MTBLS1582.html MTBLS1586.html

%.html: %.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render("'$<'", "html_document")'
