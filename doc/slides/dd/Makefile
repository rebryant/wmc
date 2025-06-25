UNIT=40pt
SCALE=0.70
PDIR = .
GEN=$(PDIR)/ddgen.py

all: 	eg-formula.tex

.SUFFIXES: .dd .idd .sdd .wdd .tex .pdf

.dd.tex:
	$(GEN) -i $< -o $@ -u $(UNIT) -s $(SCALE) -S c

.idd.tex:
	$(GEN) -i $< -o $@ -u $(UNIT) -s $(SCALE) -w -S c

.sdd.tex:
	$(GEN) -i $< -o $@ -u $(UNIT) -s 0.9 -S c

.wdd.tex:
	$(GEN) -i $< -o $@ -u $(UNIT) -w -S c

.tex.pdf:
	pdflatex $< 

eg-formula-doc.wdd: eg-formula.dd
	cp eg-formula.dd eg-formula-doc.wdd

eg-wmc-doc.wdd: eg-wmc.dd
	cp eg-wmc.dd eg-wmc-doc.wdd

decision-doc.wdd: decision.sdd
	cp decision.sdd decision-doc.wdd

product-doc.wdd: product.sdd
	cp product.sdd product-doc.wdd
