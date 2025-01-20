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

	#### Raw data
	# Variables
	$(INTERP) $(GRAB) "Data variables" mc2024_track2-random_*_*.count > all-variables-multi.csv
	# MPF Precision
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_1*.count > upos-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_2*.count > epos-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_3*.count > uposneg-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_4*.count > eposneg-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_5*.count > bposneg-mpf-prec-multi.csv
	# MPFI Est
	$(INTERP) $(GRAB) 2 "MPFI COUNT" mc2024_track2-random_*_*.count > all-mpfi-eprec.csv
	# MPFI Actual
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_3*.count > uposneg-mpfi-prec.csv
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_4*.count > eposneg-mpfi-prec.csv
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_5*.count > bposneg-mpfi-prec.csv
	# Combo times
	$(INTERP) $(GRAB) "COMBO used MPF " mc2024_track2-random_*_*.count > combo-mpf-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI " mc2024_track2-random_*_*.count > combo-mpfi-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPQ " mc2024_track2-random_*_*.count > combo-mpq-secs.csv
	# MPQ time
	$(INTERP) $(GRAB) "MPQ required " mc2024_track2-random_*_*.count > all-mpq-secs.csv
	# Minimize across runs
	$(INTERP) $(MINIMIZE) all-variables-multi.csv > all-variables.csv
	$(INTERP) $(MINIMIZE) upos-mpf-prec-multi.csv > upos-mpf-prec.csv
	$(INTERP) $(MINIMIZE) epos-mpf-prec-multi.csv > epos-mpf-prec.csv
	$(INTERP) $(MINIMIZE) uposneg-mpf-prec-multi.csv > uposneg-mpf-prec.csv
	$(INTERP) $(MINIMIZE) eposneg-mpf-prec-multi.csv > eposneg-mpf-prec.csv
	$(INTERP) $(MINIMIZE) bposneg-mpf-prec-multi.csv > bposneg-mpf-prec.csv

	#### Combinations
	# MPF Precision vs. vars
	$(INTERP) $(MERGE) -f all-variables.csv upos-mpf-prec.csv > upos-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv epos-mpf-prec.csv > epos-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv uposneg-mpf-prec.csv > uposneg-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv eposneg-mpf-prec.csv > eposneg-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv bposneg-mpf-prec.csv > bposneg-mpf+vars.csv
	# MPFI Predictive value
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv uposneg-mpfi-prec.csv > uposneg-mpfi-est+act.csv
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv eposneg-mpfi-prec.csv > eposneg-mpfi-est+act.csv
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv bposneg-mpfi-prec.csv > bposneg-mpfi-est+act.csv
	# MPQ TIMES
	$(INTERP) $(MERGE) -f combo-mpf-secs.csv all-mpq-secs.csv > combo-mpf+mpq.csv
	$(INTERP) $(MERGE) -f combo-mpfi-secs.csv all-mpq-secs.csv >  combo-mpfi+mpq.csv
	$(INTERP) $(MERGE) -f combo-mpq-secs.csv all-mpq-secs.csv > combo-mpq+mpq.csv
