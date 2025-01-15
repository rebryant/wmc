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


halfp: prod-p$(P)_half.cnf
minp: prod-p$(P)_min.cnf
subunitp: prod-p$(P)_subunit.cnf
maxp: prod-p$(P)_max.cnf
uniformp: prod-p$(P)_uniform-$(SEED).cnf
exponentialp: prod-p$(P)_exponential-$(SEED).cnf
allp: halfp minp maxp subunitp uniformp exponentialp

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



prod-p$(P).nnf:
	$(INTERP) $(GEN) -p $(P) -c prod-p$(P).nnf

prod-p$(P)_half.cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.5 > prod-p$(P)_half.cnf

prod-p$(P)_min.cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.000000001 > prod-p$(P)_min.cnf

prod-p$(P)_subunit.cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 0.999999999 > prod-p$(P)_subunit.cnf

prod-p$(P)_max.cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -p $(P) -u 1000000000 > prod-p$(P)_max.cnf

prod-p$(P)_uniform-$(SEED).cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -r u -s $(SEED) -p $(P) > prod-p$(P)_uniform-$(SEED).cnf

prod-p$(P)_exponential-$(SEED).cnf: prod-p$(P).nnf
	$(INTERP) $(GEN) -r e -s $(SEED) -p $(P) > prod-p$(P)_exponential-$(SEED).cnf


.SUFFIXES: .cnf .nnf .count

.cnf.nnf:
	$(INTERP) $(COMPILE) -f $<

.nnf.count:
	$(COUNT) $< $*.cnf

