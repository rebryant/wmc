INTERP=python3
TDIR = ~/repos/model-counting/tools
DGEN = $(TDIR)/gen-data.py
LDATA = ../benchmarks/2024-track2-public+private/logs
LPDATA = ../benchmarks/product/sweep
RDATA = data-raw
FDATA = data-formatted



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







pdata:
	rm -rf $(FDATA)
	mkdir $(FDATA)
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=darkestblue, mark=pentagon*" \
		< $(RDATA)/prodp-mpf+vars.csv > $(FDATA)/prodp-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=redorange, mark=*" \
		< $(RDATA)/upos-mpf+vars.csv > $(FDATA)/upos-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=bluegray, mark=square*" \
		< $(RDATA)/epos-mpf+vars.csv > $(FDATA)/epos-mpf+vars.tex
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=redorange, mark=*" \
		< $(RDATA)/uposneg-mpf+vars.csv > $(FDATA)/uposneg-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=bluegray, mark=square*" \
		< $(RDATA)/eposneg-mpf+vars.csv > $(FDATA)/eposneg-mpf+vars.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.3}" \
		< $(RDATA)/bposneg-mpf+vars.csv > $(FDATA)/bposneg-mpf+vars.tex
	#
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=redorange, mark=*, mark options={scale=0.8, opacity=0.7}" \
		< $(RDATA)/uposneg-mpfi-est+act.csv > $(FDATA)/uposneg-mpfi-est+act.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=bluegray, mark=square*, mark options={scale=0.8, opacity=0.7}" \
		< $(RDATA)/eposneg-mpfi-est+act.csv > $(FDATA)/eposneg-mpfi-est+act.tex
	$(INTERP) $(DGEN) -X 2 -Y 3 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.0, opacity=0.7}" \
		< $(RDATA)/bposneg-mpfi-est+act.csv > $(FDATA)/bposneg-mpfi-est+act.tex
	#
	$(INTERP) $(DGEN) -X 3 -Y 2 -O "only marks, color=redorange, mark=*, mark options={scale=0.8, opacity=0.4}" \
		< $(RDATA)/combo-mpf+mpq.csv > $(FDATA)/combo-mpf+mpq.tex
	$(INTERP) $(DGEN) -X 3 -Y 2 -O "only marks, color=bluegray, mark=square*, mark options={scale=0.8, opacity=0.4}" \
		< $(RDATA)/combo-mpfi+mpq.csv > $(FDATA)/combo-mpfi+mpq.tex
	$(INTERP) $(DGEN) -X 3 -Y 2 -O "only marks, color=darkgreen, mark=diamond*, mark options={scale=1.0, opacity=0.4}" \
		< $(RDATA)/combo-mpq+mpq.csv > $(FDATA)/combo-mpq+mpq.tex

