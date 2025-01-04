INTERP = python3
COMPILE = ../../tools/compile.py
COUNT = ../../bin/nnfcount
GEN = ../generators/max_precision.py
N=10
P=1


precision-$(N).nnf:
	$(INTERP) $(GEN) -n $(N) -r precision-$(N)

precision-10pwr$(P).nnf:
	$(INTERP) $(GEN) -N $(P) -r precision-10pwr$(P)


.SUFFIXES: .nnf .scount

.nnf.scount:
	$(COUNT) -s $< $*.cnf

