all: MTBLS1582.html MTBLS1586.html Coverage.html

%.html: %.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render("'$<'", "html_document")'
