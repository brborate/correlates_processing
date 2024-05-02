## data_processed         : create processed data from raw data
data_processed: check_raw_data risk_analysis make_clean_data check_clean_data deploy_processed_dataset

check_raw_data:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_boost id27hpv))
else 
	Rscript data_clean/make_raw_dat_check.R
endif


risk_analysis:  
ifeq ($(TRIAL),$(filter $(TRIAL), id27hpv covail nvx_uk302 prevent19_stage2 azd1222_stage2))
else
	$(MAKE) -k -C riskscore_baseline all
endif


TARGET_FILE := data_clean/make_dat_$(TRIAL).R

make_clean_data: 
ifneq ($(wildcard $(TARGET_FILE)),)
	Rscript $(TARGET_FILE)
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_partA_VL))
	Rscript data_clean/RunhotdeckMI_janssen_partA_VL.R
	Rscript data_clean/make_dat_proc.R
else ifeq ($(TRIAL),$(filter $(TRIAL), vat08_combined))
	Rscript data_clean/RunhotdeckMI_sanofi_bothtrials_PartA.R
	Rscript data_clean/make_dat_vat08_combined.R
else 
	Rscript data_clean/make_dat_proc.R
endif


check_clean_data: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_boost id27hpv covail))
else 
	Rscript data_clean/make_clean_dat_check.R
endif	


## risk_report            : builds the CoVPN baseline risk score report
risk_report: data_processed
	    bash ./_build.sh riskscore

## risk_report_clean      : cleans and builds the CoVPN baseline risk score report
risk_report_clean: risk_clean data_processed
	    bash ./_build.sh riskscore

## risk_clean      : wipes out TRIAL directory within riskscore_baseline/output
risk_clean:
	    rm -rf riskscore_baseline/output/$(TRIAL)/*



## deploy risk score processed dataset : Deploy risk score dataset after checking
deploy_processed_dataset:
	    Rscript data_clean/deploy_risk_score_dataset.R





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

.PHONY: style help data_processed
