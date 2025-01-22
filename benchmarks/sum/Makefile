INTERP = python3
HOME=/Users/bryant/repos/wmc
COMPILE = $(HOME)/tools/compile.py
COUNT = $(HOME)/bin/nnfcount
GEN = $(HOME)/benchmarks/generators/gen_sum.py
N=10
P=1
SEED=123456

half: half-$(N).cnf
min: min-$(N).cnf
max: max-$(N).cnf
subunit: subunit-$(N).cnf
uniform: uniform-$(N)-$(SEED).cnf
exponential: exponential-$(N)-$(SEED).cnf


sump: sump-$(P).cnf

half-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.5 > half-$(N).cnf

min-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.000000001 > min-$(N).cnf

max-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 1000000000 > max-$(N).cnf

subunit-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.999999999 > subunit-$(N).cnf

uniform-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r u -s $(SEED) -n $(N) > uniform-$(N)-$(SEED).cnf

exponential-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r e -s $(SEED) -n $(N) > exponential-$(N)-$(SEED).cnf



sump-$(P).cnf:
	$(INTERP) $(GEN) -p $(P) -u 1.0 > sump-$(P).cnf

.SUFFIXES: .cnf .nnf .count

.cnf.nnf:
	$(INTERP) $(COMPILE) -f $<

.nnf.count:
	$(COUNT) $< $*.cnf

