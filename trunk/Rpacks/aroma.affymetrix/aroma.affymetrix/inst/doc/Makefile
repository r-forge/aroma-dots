foo = aroma.affymetrix-DevelopersGuide

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Includes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
BMCINCLUDES = bmc_article.cls bmc_article.bst
BIBINCLUDES=hb-at-maths.lth.se.bib 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Application settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gs=gs
gsview=gsview32

#gs="C:/gs/gs8.51/bin/gswin32c.exe"
#gsview="C:/Program Files/GSView/gsview/gsview32.exe"



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Final
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$(foo).ps: $(foo).dvi
	dvips -Ppdf -ta4 $(foo).dvi

$(foo).pdf: $(foo).ps
	ps2pdf -sPAPERSIZE=a4 $(foo).ps

$(foo).dvi: $(foo).tex $(foo).bbl ${BIBINCLUDES}
	latex $(foo).tex; 
	latex $(foo).tex

$(foo).bbl: $(foo).tex $(foo).bib 
	latex $(foo).tex;
	bibtex $(foo);
	latex $(foo).tex;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bibliography
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
$(foo).bib: ${BIBINCLUDES}
	grep -v "url[ ]*=" $(BIBINCLUDES) > $(foo).bib


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spell check submit version
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
spell: $(foo).tex
	export TERM=cygwin
	aspell --mode=tex --check --lang=EN $(foo).tex

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Shortcuts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
bib: $(foo).bib
dvi: $(foo).dvi
ps: $(foo).ps
pdf: $(foo).pdf

src: $(foo)/
src.zip: $(foo).zip
src.tar: $(foo).tar
src.tar.gz: $(foo).tar.gz 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean up
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
clean:
	rm -f tmp.ps;
	rm -f $(foo)*.aux;
	rm -f $(foo)*.blg;
	rm -f $(foo)*.bbl;
	rm -f $(foo)*.dvi;
	rm -f $(foo)*.fff;
	rm -f $(foo)*.log;
	rm -f $(foo)*.nav;
	rm -f $(foo)*.out;
	rm -f $(foo)*.snm;
	rm -f $(foo)*.toc;
	rm -f $(foo)*.ttt;
	rm -f $(foo)*.vrb;
