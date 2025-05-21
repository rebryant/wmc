INTERP=python3
TDIR = ~/repos/model-counting/tools
DGEN = $(TDIR)/gen-data.py
LDATA = ../benchmarks/2024-track2-public+private/logs
LPDATA = ../benchmarks/optimized-product
RDATA = data-raw
FDATA = data-formatted
TABULATE = tabulate_combo.py
OPTIMIZE = optimize_performance.py


rdata:
	rm -rf $(RDATA)
	mkdir $(RDATA)
	cp -p $(LPDATA)/optimized-product-1000.csv $(RDATA)
	cp -p $(LDATA)/upos-mpf+vars.csv $(RDATA)
	cp -p $(LDATA)/epos-mpf+vars.csv $(RDATA)
	cp -p $(LDATA)/uposneg-mpf+vars.csv $(RDATA)
	cp -p $(LDATA)/eposneg-mpf+vars.csv $(RDATA)
	cp -p $(LDATA)/bposneg-mpf+vars.csv $(RDATA)
	cp -p $(LDATA)/uposneg-mpfi-est+act.csv $(RDATA)
	cp -p $(LDATA)/eposneg-mpfi-est+act.csv $(RDATA)
	cp -p $(LDATA)/bposneg-mpfi-est+act.csv $(RDATA)
	cp -p $(LDATA)/combo-mpf+mpq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpfi+mpq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpfi2+mpq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpq+mpq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpf+nompq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpfi+nompq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpfi2+nompq.csv $(RDATA)
	cp -p $(LDATA)/tabulate-all.csv $(RDATA)
	$(INTERP) $(TABULATE) $(RDATA) > method_table.tex
	cp -p $(LDATA)/negpos-tabulate.csv $(RDATA)
	cp -p $(LDATA)/nonneg-tabulate.csv $(RDATA)



pdata:
	rm -rf $(FDATA)
	mkdir $(FDATA)
	#
	$(INTERP) $(DGEN) -X 1 -Y 3 -y 50 -O "only marks, mark=pentagon*, fill=blue!10, color=darkestblue, mark options={scale=0.8, opacity=0.7}" \
		< $(RDATA)/optimized-product-1000.csv > $(FDATA)/optimized-product.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, color=redorange, mark=*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/upos-mpf+vars.csv > $(FDATA)/upos-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/epos-mpf+vars.csv > $(FDATA)/epos-mpf+vars.tex
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, color=redorange, mark=*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/uposneg-mpf+vars.csv > $(FDATA)/uposneg-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/eposneg-mpf+vars.csv > $(FDATA)/eposneg-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.1, opacity=0.7}" \
		< $(RDATA)/bposneg-mpf+vars.csv > $(FDATA)/bposneg-mpf+vars.tex
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -x 50 -y 50 -O "only marks, color=redorange, mark=*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/uposneg-mpfi-est+act.csv > $(FDATA)/uposneg-mpfi-est+act.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -x 50 -y 50  -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/eposneg-mpfi-est+act.csv > $(FDATA)/eposneg-mpfi-est+act.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -x 50 -y 50 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.1, opacity=0.7}" \
		< $(RDATA)/bposneg-mpfi-est+act.csv > $(FDATA)/bposneg-mpfi-est+act.tex
	#
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=redorange, mark=*, mark options={scale=0.6, opacity=0.7}" \
		< $(RDATA)/combo-mpf+mpq.csv > $(FDATA)/combo-mpf+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/combo-mpfi+mpq.csv > $(FDATA)/combo-mpfi+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 4 -Y 5 -O "only marks, color=midblue, mark=diamond*, mark options={scale=0.9, opacity=0.7}" \
		< $(RDATA)/combo-mpfi2+mpq.csv > $(FDATA)/combo-mpfi2+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 4 -Y 5 -O "only marks, color=darkgreen, mark=pentagon*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/combo-mpq+mpq.csv > $(FDATA)/combo-mpq+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=redorange, mark=*, mark options={scale=0.6, opacity=0.7}" \
		< $(RDATA)/combo-mpf+nompq.csv > $(FDATA)/combo-mpf+nompq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/combo-mpfi+nompq.csv > $(FDATA)/combo-mpfi+nompq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 4 -Y 5 -O "only marks, color=midblue, mark=diamond*, mark options={scale=0.9, opacity=0.7}" \
		< $(RDATA)/combo-mpfi2+nompq.csv > $(FDATA)/combo-mpfi2+nompq.tex
	$(INTERP) $(OPTIMIZE) -l -d $(RDATA) -m c > $(FDATA)/tabulate-count-limited.tex
	$(INTERP) $(OPTIMIZE) -l -d $(RDATA) -m e > $(FDATA)/tabulate-effort-limited.tex
