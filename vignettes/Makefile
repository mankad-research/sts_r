# Compress the vignette with ghostscript when it is available (gs, or
# gswin64c/gswin32c on Windows); otherwise copy it uncompressed.
GS := $(firstword $(foreach g,gs gswin64c gswin32c,$(shell command -v $(g) 2>/dev/null)))

all: ../inst/doc/vignette.pdf

../inst/doc/vignette.pdf: usingSTS.pdf
	mkdir -p ../inst/doc
ifneq ($(GS),)
	"$(GS)" -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
		-dPDFSETTINGS=/ebook \
		-sOutputFile=../inst/doc/vignette.pdf usingSTS.pdf
else
	cp usingSTS.pdf ../inst/doc/vignette.pdf
endif
