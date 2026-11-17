PAPER = ms
SUBDIR = submit
ARXIVDIR = arxiv

.PHONY: all clean submit arxiv

all:
	lualatex $(PAPER)
	bibtex $(PAPER)
	lualatex $(PAPER)
	lualatex $(PAPER)

# Flat (no subdirectories) bundle for the AAS upload system: build the paper,
# copy its sources into $(SUBDIR)/, zip them, then build the flat copy in
# $(SUBDIR)-check/ to confirm it compiles without figures/ and tables/.
submit: all
	rm -rf $(SUBDIR) $(SUBDIR)-check $(SUBDIR).zip
	mkdir $(SUBDIR)
	cp $(PAPER).tex $(PAPER).bbl refs.bib aastex702.cls aasjournalv7.1.bst \
	   figures/*.pdf tables/*.tex $(SUBDIR)/
	cd $(SUBDIR) && zip -q ../$(SUBDIR).zip *
	cp -R $(SUBDIR) $(SUBDIR)-check
	cd $(SUBDIR)-check && lualatex -interaction=nonstopmode -halt-on-error $(PAPER) > /dev/null \
	   && bibtex $(PAPER) > /dev/null \
	   && lualatex -interaction=nonstopmode -halt-on-error $(PAPER) > /dev/null \
	   && lualatex -interaction=nonstopmode -halt-on-error $(PAPER) > /dev/null
	@echo "Wrote $(SUBDIR).zip ($$(ls $(SUBDIR) | wc -l | tr -d ' ') files); flat build OK: $(SUBDIR)-check/$(PAPER).pdf"

# Flat arXiv bundle: arXiv compiles with pdfLaTeX and does not run BibTeX, so
# ship the .bbl (not refs.bib or the .bst), then test-build the flat copy in
# $(ARXIVDIR)-check/ with pdflatex alone, the way arXiv will.
arxiv: all
	rm -rf $(ARXIVDIR) $(ARXIVDIR)-check $(ARXIVDIR).tar.gz
	mkdir $(ARXIVDIR)
	cp $(PAPER).bbl aastex702.cls figures/*.pdf tables/*.tex $(ARXIVDIR)/
	sed -e '1s/,linenumbers//' -e '1s/linenumbers,//' $(PAPER).tex > $(ARXIVDIR)/$(PAPER).tex
	@! grep -q 'linenumbers' $(ARXIVDIR)/$(PAPER).tex || { echo "linenumbers still in $(ARXIVDIR)/$(PAPER).tex"; exit 1; }
	cd $(ARXIVDIR) && tar -czf ../$(ARXIVDIR).tar.gz *
	cp -R $(ARXIVDIR) $(ARXIVDIR)-check
	cd $(ARXIVDIR)-check && for i in 1 2 3; do \
	   pdflatex -interaction=nonstopmode -halt-on-error $(PAPER) > /dev/null || exit 1; done
	@grep -E 'undefined|Rerun' $(ARXIVDIR)-check/$(PAPER).log || true
	@echo "Wrote $(ARXIVDIR).tar.gz ($$(ls $(ARXIVDIR) | wc -l | tr -d ' ') files); pdflatex build OK: $(ARXIVDIR)-check/$(PAPER).pdf"

clean:
	rm -f $(shell ls $(PAPER).* | grep -Ev '\.(tex|bib|pdf)$$')
	rm -rf $(SUBDIR) $(SUBDIR)-check $(SUBDIR).zip
	rm -rf $(ARXIVDIR) $(ARXIVDIR)-check $(ARXIVDIR).tar.gz
