INTERP=python3
TDIR = ~/repos/model-counting/tools
DGEN = $(TDIR)/gen-data.py
LDATA = ../benchmarks/2024-track2-public+private/logs
LPDATA = ../benchmarks/product/sweep
RDATA = data-raw
FDATA = data-formatted
TABULATE = tabulate_combo.py


rdata:
	rm -rf $(RDATA)
	mkdir $(RDATA)
	cp -p $(LPDATA)/prodp-mpf+vars.csv $(RDATA)
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
	cp -p $(LDATA)/combo-mpq+mpq.csv $(RDATA)
	cp -p $(LDATA)/combo-mpf-only-secs.csv $(RDATA)
	cp -p $(LDATA)/combo-mpfi-only-secs.csv $(RDATA)
	$(INTERP) $(TABULATE) $(RDATA) > method_table.tex



pdata:
	rm -rf $(FDATA)
	mkdir $(FDATA)
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -y 50 -O "only marks, mark=pentagon*, fill=blue!10, color=darkestblue, mark options={scale=0.8, opacity=0.7}" \
		< $(RDATA)/prodp-mpf+vars.csv > $(FDATA)/prodp-mpf+vars.tex
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
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=redorange, mark=*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/combo-mpf+mpq.csv > $(FDATA)/combo-mpf+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=bluegray, mark=triangle*, mark options={scale=0.7, opacity=0.7}" \
		< $(RDATA)/combo-mpfi+mpq.csv > $(FDATA)/combo-mpfi+mpq.tex
	$(INTERP) $(DGEN) -L 0.001 -X 3 -Y 2 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.1, opacity=0.7}" \
		< $(RDATA)/combo-mpq+mpq.csv > $(FDATA)/combo-mpq+mpq.tex

