TDIR = ~/repos/model-counting/tools
INTERP=python3
GRAB = $(TDIR)/grab_data.py
MERGE = $(TDIR)/merge_csv.py
EXTRACT = $(TDIR)/extract_columns.py
LTDIR = ~/repos/wmc/tools
MINIMIZE = $(LTDIR)/combine_data.py


clean:
	rm -f *~
	rm -f *.csv

data:

	$(INTERP) $(GRAB) "Data variables" prodp-*_*.count > all-variables-multi.csv
	# MPF precision
	$(INTERP) $(GRAB) 2 "MPF COUNT" prodp-*_*.count > prodp-mpf-multi.csv
	# Minimize across runs
	$(INTERP) $(MINIMIZE) all-variables-multi.csv > all-variables.csv
	$(INTERP) $(MINIMIZE) prodp-mpf-multi.csv > prodp-mpf.csv
	$(INTERP) $(MERGE) -f all-variables.csv prodp-mpf.csv > prodp-mpf+vars.csv

