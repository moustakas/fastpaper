PAPER = ms
SUBDIR = submit

.PHONY: all clean submit

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

clean:
	rm -f $(shell ls $(PAPER).* | grep -Ev '\.(tex|bib|pdf)$$')
	rm -rf $(SUBDIR) $(SUBDIR)-check $(SUBDIR).zip
