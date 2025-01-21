INTERP = python3
HOME=/Users/bryant/repos/wmc
COMPILE = $(HOME)/tools/compile.py
COUNT = $(HOME)/bin/nnfcount
GEN = $(HOME)/benchmarks/generators/gen_product.py
N=10
P=1
SEED=123456

half: half-$(N).cnf
min: min-$(N).cnf
min: max-$(N).cnf
uniform: uniform-$(N)-$(SEED).cnf
exponential: exponential-$(N)-$(SEED).cnf


prodp: prodp-$(P).cnf
halfp: prodp-$(P)_half.cnf
minp: prodp-$(P)_min.cnf
subunitp: prodp-$(P)_subunit.cnf
hugep: prodp-$(P)_huge.cnf
maxp: prodp-$(P)_max.cnf
uniformp: prodp-$(P)_uniform-$(SEED).cnf
exponentialp: prodp-$(P)_exponential-$(SEED).cnf
allp: prodp halfp hugep minp maxp subunitp uniformp exponentialp

half-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.5 > half-$(N).cnf

min-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.000000001 > min-$(N).cnf

max-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 1000000000 > max-$(N).cnf

subunit-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.999999999 > subunit-$(N).cnf

pessimistic-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) > pessimistic-$(N).cnf

uniform-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r u -s $(SEED) -n $(N) > uniform-$(N)-$(SEED).cnf

exponential-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r e -s $(SEED) -n $(N) > exponential-$(N)-$(SEED).cnf


prodp-$(P).nnf:
	$(INTERP) $(GEN) -p $(P) -c prodp-$(P).nnf

prodp-$(P).cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 1.0 > prodp-$(P).cnf

prodp-$(P)_half.cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.5 > prodp-$(P)_half.cnf

prodp-$(P)_min.cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.000000001 > prodp-$(P)_min.cnf

prodp-$(P)_subunit.cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.999999999 > prodp-$(P)_subunit.cnf

prodp-$(P)_max.cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 1000000000 > prodp-$(P)_max.cnf

prodp-$(P)_huge.cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u  999999999.999999999 > prodp-$(P)_huge.cnf

prodp-$(P)_uniform-$(SEED).cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -r u -s $(SEED) -p $(P) > prodp-$(P)_uniform-$(SEED).cnf

prodp-$(P)_exponential-$(SEED).cnf: prodp-$(P).nnf
	$(INTERP) $(GEN) -r e -s $(SEED) -p $(P) > prodp-$(P)_exponential-$(SEED).cnf


.SUFFIXES: .cnf .nnf .count

.cnf.nnf:
	$(INTERP) $(COMPILE) -f $<

.nnf.count:
	$(COUNT) $< $*.cnf

