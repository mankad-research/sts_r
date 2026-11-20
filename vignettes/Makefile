# Compress the vignette with ghostscript when it is available (gs, or
# gswin64c/gswin32c on Windows); otherwise copy it uncompressed. The lookup
# runs in the shell so that paths with spaces ("C:/Program Files/...") work.

all: ../inst/doc/vignette.pdf

../inst/doc/vignette.pdf: usingSTS.pdf
	mkdir -p ../inst/doc
	GS="$$(command -v gs || command -v gswin64c || command -v gswin32c || true)"; \
	if [ -n "$$GS" ]; then \
		"$$GS" -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
			-dPDFSETTINGS=/ebook \
			-sOutputFile=../inst/doc/vignette.pdf usingSTS.pdf; \
	else \
		cp usingSTS.pdf ../inst/doc/vignette.pdf; \
	fi
