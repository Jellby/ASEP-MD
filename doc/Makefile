ifeq ($(LANG_),English)
  XSLFLAG = -stringparam lang en
else
  ifeq ($(LANG_),Espanol)
    XSLFLAG = -stringparam lang es
  else
  endif
endif

XSLT   = xsltproc $(XSLFLAG)
RMF    = rm -f

FUNCS = functions.xsl
STYLE = style.sty
XML = method.xml input.xml output.xml

TRASH = $(XML:.xml=.aux) $(XML:.xml=.log) x.log

#==============================================================================

.PHONY: all xhtml pdf clean

.SECONDARY: $(XML:.xml=.tex)

all: txt xhtml pdf

txt: $(XML:.xml=.txt)

xhtml: $(XML:.xml=.xhtml)

pdf: $(XML:.xml=.pdf)

%.txt: %.xml %_txt.xsl txt.xsl $(FUNCS)
	$(XSLT) $*_txt.xsl $*.xml > $@

%.xhtml: %.xml %_xhtml.xsl xhtml.xsl $(FUNCS)
	$(XSLT) $*_xhtml.xsl $*.xml > $@

%.tex: %.xml %_tex.xsl tex.xsl $(FUNCS)
	$(XSLT) $*_tex.xsl $*.xml > $@

%.pdf: %.tex $(STYLE)
	pdflatex $*.tex
	pdflatex $*.tex
	$(RMF) $*.aux $*.log x.log

clean:
	$(RMF) $(XML:.xml=.txt) $(XML:.xml=.xhtml) $(XML:.xml=.tex) $(XML:.xml=.pdf) $(TRASH)
