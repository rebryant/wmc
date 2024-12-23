INTERP = python3
COMPILE = ../../tools/compile.py
COUNT = ../../bin/nnfcount
GEN = ../generators/gen_product.py
N=10
SEED=123456

half: half-$(N).cnf
min: min-$(N).cnf
min: max-$(N).cnf
pessimistic: pessimistic-$(N).cnf
uniform: uniform-$(N)-$(SEED).cnf
exponential: exponential-$(N)-$(SEED).cnf

half-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.5 > half-$(N).cnf

min-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.000000001 > min-$(N).cnf

max-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) -u 0.999999999 > max-$(N).cnf

pessimistic-$(N).cnf:
	$(INTERP) $(GEN) -n $(N) > pessimistic-$(N).cnf

uniform-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r u -s $(SEED) -n $(N) > uniform-$(N)-$(SEED).cnf

exponential-$(N)-$(SEED).cnf:
	$(INTERP) $(GEN) -r e -s $(SEED) -n $(N) > exponential-$(N)-$(SEED).cnf


.SUFFIXES: .cnf .nnf .count

.cnf.nnf:
	$(INTERP) $(COMPILE) -f $<

.nnf.count:
	$(COUNT) $< $*.cnf

