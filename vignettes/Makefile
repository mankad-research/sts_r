all: ../inst/doc/vignette.pdf

../inst/doc/vignette.pdf: usingSTS.pdf
	mkdir -p ../inst/doc
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
		-dPDFSETTINGS=/ebook \
		-sOutputFile=../inst/doc/vignette.pdf usingSTS.pdf