TDIR = ~/repos/model-counting/tools
INTERP=python3
GRAB = $(TDIR)/grab_data.py
MERGE = $(TDIR)/merge_csv.py
EXTRACT = $(TDIR)/extract_columns.py
LTDIR = ~/repos/wmc/tools
MINIMIZE = $(LTDIR)/combine_data.py
PREC256-LIMITED = prec256-postmpfi
PREC256-ALL = prec256-all
ORIGINAL = original

clean:
	rm -f *~
	rm -f *.csv

data:

	#### Raw data
	# Variables
	$(INTERP) $(GRAB) "Data variables" mc2024_track2-random_*_*.count > all-variables-multi.csv
	# Operations
	$(INTERP) $(GRAB) "Binops  TOTAL" mc2024_track2-random_*_*.count > all-binops-multi.csv
	# MPF Precision
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_1*.count > upos-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_2*.count > epos-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_3*.count > uposneg-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_4*.count > eposneg-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" mc2024_track2-random_*_5*.count > bposneg-mpf-prec-multi.csv
	$(INTERP) $(GRAB) 2 "MPF COUNT" $(ORIGINAL)/mc2024_track2-random_*.count > original-mpf-prec.csv
	# MPFI Est
	$(INTERP) $(GRAB) 2 "MPFI COUNT" mc2024_track2-random_*_*.count > all-mpfi-eprec.csv
	# MPFI Actual
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_3*.count > uposneg-mpfi-prec.csv
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_4*.count > eposneg-mpfi-prec.csv
	$(INTERP) $(GRAB) 3 "MPFI COUNT" mc2024_track2-random_*_5*.count > bposneg-mpfi-prec.csv
	# Combo times
	$(INTERP) $(GRAB) "COMBO used MPF " mc2024_track2-random_*_*.count > combo-mpf-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI " mc2024_track2-random_*_*.count > combo-mpfi-secs.csv
	$(INTERP) $(GRAB) "MPFI gave only" mc2024_track2-random_*_*.count > combo-mpfi-fail-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPQ " mc2024_track2-random_*_*.count > combo-mpq-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPQ_ABORT" mc2024_track2-random_*_*.count > combo-mpfi-abort-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPF_ONLY" mc2024_track2-random_*_*.count > combo-mpf-only-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI_ONLY" mc2024_track2-random_*_*.count > combo-mpfi-only-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI" $(PREC256-LIMITED)/mc2024_track2-random_*_*.count > combo-mpfi2-postmpfi-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI" $(PREC256-ALL)/mc2024_track2-random_*_*.count > combo-mpfi2-secs.csv
	$(INTERP) $(GRAB) "MPFI gave only" $(PREC256-LIMITED)/mc2024_track2-random_*_*.count > combo-mpfi2-fail-secs.csv
	$(INTERP) $(GRAB) "COMBO used MPFI_ONLY" $(PREC256-LIMITED)/mc2024_track2-random_{023,178}_*.count > combo-mpfi2-only-secs.csv
	# Times
	$(INTERP) $(GRAB) "MPF required " mc2024_track2-random_*_*.count > normal-mpf-secs.csv
	$(INTERP) $(GRAB) "MPFI required " mc2024_track2-random_*_*.count > normal-mpfi-secs.csv
	$(INTERP) $(GRAB) "MPQ required " mc2024_track2-random_*_*.count > all-mpq-secs.csv
	cat normal-mpf-secs.csv combo-mpf-only-secs.csv | sort > all-mpf-secs.csv
	cat normal-mpfi-secs.csv combo-mpfi-only-secs.csv combo-mpfi-abort-secs.csv | sort > all-mpfi-secs.csv
	cat combo-mpfi2-secs.csv combo-mpfi2-fail-secs.csv | sort > all-mpfi2-secs.csv
	# Minimize across runs
	$(INTERP) $(MINIMIZE) all-variables-multi.csv > all-variables.csv
	$(INTERP) $(MINIMIZE) all-binops-multi.csv > all-binops.csv
	$(INTERP) $(MINIMIZE) upos-mpf-prec-multi.csv > upos-mpf-prec.csv
	$(INTERP) $(MINIMIZE) epos-mpf-prec-multi.csv > epos-mpf-prec.csv
	$(INTERP) $(MINIMIZE) uposneg-mpf-prec-multi.csv > uposneg-mpf-prec.csv
	$(INTERP) $(MINIMIZE) eposneg-mpf-prec-multi.csv > eposneg-mpf-prec.csv
	$(INTERP) $(MINIMIZE) bposneg-mpf-prec-multi.csv > bposneg-mpf-prec.csv

	#### Combinations
	# d-DNNF information
	$(INTERP) $(MERGE) all-variables.csv all-binops.csv > vars+binops.csv
	# MPF Precision vs. vars
	$(INTERP) $(MERGE) -f all-variables.csv upos-mpf-prec.csv > upos-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv epos-mpf-prec.csv > epos-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv uposneg-mpf-prec.csv > uposneg-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv eposneg-mpf-prec.csv > eposneg-mpf+vars.csv
	$(INTERP) $(MERGE) -f all-variables.csv bposneg-mpf-prec.csv > bposneg-mpf+vars.csv
	$(INTERP) $(MERGE) -f -t 2:49:49 all-variables.csv original-mpf-prec.csv > original-mpf+vars.csv
	# MPFI Predictive value
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv uposneg-mpfi-prec.csv > uposneg-mpfi-est+act.csv
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv eposneg-mpfi-prec.csv > eposneg-mpfi-est+act.csv
	$(INTERP) $(MERGE) -f all-mpfi-eprec.csv bposneg-mpfi-prec.csv > bposneg-mpfi-est+act.csv
	# MPQ TIMES
	$(INTERP) $(MERGE) -f combo-mpf-secs.csv all-mpq-secs.csv > combo-mpf+mpq.csv
	$(INTERP) $(MERGE) -f combo-mpfi-secs.csv all-mpq-secs.csv >  combo-mpfi+mpq.csv
	$(INTERP) $(MERGE) -s 1:2 -f combo-mpfi-fail-secs.csv combo-mpfi2-secs.csv all-mpq-secs.csv > combo-mpfi2+mpq.csv
	$(INTERP) $(MERGE) -s 1:2:3 -f combo-mpfi-fail-secs.csv combo-mpfi2-fail-secs.csv all-mpq-secs.csv > combo-mpq+mpq.csv
	$(INTERP) $(MERGE) combo-mpf-only-secs.csv 13000 > combo-mpf+nompq.csv
	$(INTERP) $(MERGE) combo-mpfi-only-secs.csv 13000 > combo-mpfi+nompq.csv
	$(INTERP) $(MERGE) -f -s 1:2 combo-mpfi2-only-secs.csv combo-mpfi-fail-secs.csv 13000 > combo-mpfi2+nompq.csv
	# Modalities
	cat combo-mpf-secs.csv combo-mpf-only-secs.csv | sort > used-mpf-secs.csv
	cat combo-mpfi-secs.csv combo-mpfi-only-secs.csv | sort > used-mpfi-secs.csv
	$(INTERP) $(MERGE) -t -1:5000:1 used-mpf-secs.csv > modality-mpf.csv
	$(INTERP) $(MERGE) -t -1:5000:2 used-mpfi-secs.csv > modality-mpfi.csv
	$(INTERP) $(MERGE) -t -1:5000:3 combo-mpfi2-postmpfi-secs.csv > modality-mpfi2.csv
	$(INTERP) $(MERGE) -t -1:5000:4 combo-mpfi2-fail-secs.csv > modality-mpq.csv
	cat modality-mpf.csv modality-mpfi.csv modality-mpfi2.csv modality-mpq.csv | sort | $(INTERP) $(EXTRACT) 1:3 > modality-all.csv
	$(INTERP) $(MERGE) -l bench,mode,mpf,mpfi,mpfi2,mpq modality-all.csv all-mpf-secs.csv all-mpfi-secs.csv all-mpfi2-secs.csv all-mpq-secs.csv > tabulate-all.csv
	# Nonegative using MPF
	$(INTERP) $(GRAB) "Data variables" mc2024_track2-random_*_{1,2}*.count > nonneg-variables.csv	
#	$(INTERP) $(GRAB) 2 "WEIGHTED DBL COUNT" mc2024_track2-random_*_{1,2}*.count > nonneg-double-precision.csv	
#	$(INTERP) $(GRAB) "DBL required" mc2024_track2-random_*_{1,2}*.count > nonneg-double-seconds.csv
	$(INTERP) $(GRAB) "COMBO used ERD_ONLY" prec-erd/mc2024_track2-random_*_{1,2}*.count > nonneg-erd-seconds.csv
	$(INTERP) $(GRAB) "COMBO used MPF_ONLY" prec64-all/mc2024_track2-random_*_{1,2}*.count > nonneg-mpf64-seconds.csv
	$(INTERP) $(GRAB) "COMBO used MPF" mc2024_track2-random_*_{1,2}*.count > nonneg-mpf128-seconds.csv
	$(INTERP) $(GRAB) "COMBO used MPF_ONLY" prec256-all/mc2024_track2-random_*_{1,2}*.count > nonneg-mpf256-seconds.csv
	$(INTERP) $(GRAB) "MPQ required " mc2024_track2-random_*_{1,2}*.count > nonneg-mpq-seconds.csv
	$(INTERP) $(MERGE) -l bench,var,erd-sec,mpf64-sec,mpf128-sec,mpf256-sec,mpq-sec \
		nonneg-variables.csv nonneg-erd-seconds.csv nonneg-mpf64-seconds.csv  nonneg-mpf128-seconds.csv nonneg-mpf256-seconds.csv nonneg-mpq-seconds.csv > nonneg-tabulate.csv
	# Negative/positive using MPFI
	$(INTERP) $(GRAB) "Data variables" mc2024_track2-random_*_{3,4,5}*.count > negpos-variables.csv
	$(INTERP) $(GRAB) 2 "COMBO COUNT" prec64-all/mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi064-precision.csv
	$(INTERP) $(GRAB) "COMBO used" prec64-all/mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi064-seconds.csv
	$(INTERP) $(GRAB) 2 "WEIGHTED MPFI COUNT" mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi128-precision-normal.csv
	$(INTERP) $(GRAB) "MPFI required" mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi128-seconds-normal.csv
	$(INTERP) $(GRAB) 2 "COMBO COUNT" mc2024_track2-random_{023,178}_{3,4,5}*.count > negpos-mpfi128-precision-hard.csv
	$(INTERP) $(GRAB) "COMBO used" mc2024_track2-random_{023,178}_{3,4,5}*.count > negpos-mpfi128-seconds-hard.csv
	cat negpos-mpfi128-precision-normal.csv negpos-mpfi128-precision-hard.csv | sort > negpos-mpfi128-precision.csv
	cat negpos-mpfi128-seconds-normal.csv negpos-mpfi128-seconds-hard.csv | sort > negpos-mpfi128-seconds.csv
	$(INTERP) $(GRAB) 2 "COMBO COUNT" prec256-all/mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi256-precision.csv
	$(INTERP) $(GRAB) "COMBO used" prec256-all/mc2024_track2-random_*_{3,4,5}*.count > negpos-mpfi256-seconds.csv
	$(INTERP) $(GRAB) "MPQ required " mc2024_track2-random_*_{3,4,5}*.count > negpos-mpq-seconds.csv
	$(INTERP) $(MERGE) -l bench,var,mpfi64-prec,mpfi64-sec,mpfi128-prec,mpfi128-sec,mpfi256-prec,mpfi256-sec,mpq-sec \
		negpos-variables.csv negpos-mpfi064-precision.csv negpos-mpfi064-seconds.csv negpos-mpfi128-precision.csv negpos-mpfi128-seconds.csv negpos-mpfi256-precision.csv negpos-mpfi256-seconds.csv negpos-mpq-seconds.csv > negpos-tabulate.csv
