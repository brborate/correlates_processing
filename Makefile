## immuno_analysis        : builds immunogenicity exploratory analyses
immuno_analysis: 
	$(MAKE) -k -C immuno_tabular all
	$(MAKE) -k -C immuno_graphical all

## immuno_report          : builds the CoVPN immunogenicity report
immuno_report: immuno_analysis
	bash ./_build.sh immuno

## risk_analysis          : builds Baseline Risk Score analysis
risk_analysis: data_processed
	$(MAKE) -k -C riskscore_baseline all

## risk_report            : builds the CoVPN baseline risk score report
risk_report: risk_analysis
	bash ./_build.sh riskscore

## data_processed         : create processed data from raw data
data_processed: check_raw_data make_clean_data check_clean_data

check_raw_data:
	Rscript data_clean/make_raw_dat_check.R
make_clean_data: check_raw_data
	Rscript data_clean/make_dat_proc.R
check_clean_data: make_clean_data
	Rscript data_clean/make_clean_dat_check.R
## help_checks            : see a list of checks that are run on the data during cleaning
help_tests: data_clean/make_clean_dat_check.R data_clean/make_raw_dat_check.R
	@echo "\nTests on the raw data: \n"
	@sed -n 's/^##//p' data_clean/make_raw_dat_check.R
	@echo "\nTests on the clean data: \n"
	@sed -n 's/^##//p' data_clean/make_clean_dat_check.R
	@echo "\n"

## style                  : re-styles the codebase for consistent formatting
style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: style help immuno_analysis \
  immuno_report cor_report cor_analysis data_processed
