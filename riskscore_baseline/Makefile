## all           : clean saved results and get results to populate report
all: get_riskscores

## clean         : delete all SL results
clean:
	rm -rf output/$(TRIAL)/*

## get_riskscores : append previous risk scores if input data has not changed OR generate new risk scores!
get_riskscores:  
ifeq ($(TRIAL), $(filter $(TRIAL), vat08m vat08b))
	Rscript code/get_riskscores.R bseroneg
	Rscript code/get_riskscores.R bseropos
	Rscript code/stack_bseroneg_bseropos.R
ifeq ($(TRIAL), $(filter $(TRIAL), moderna_boost))
	Rscript code/get_riskscores_for_COVEBoost.R
else
	Rscript code/get_riskscores.R 
endif

## run_cvsl_riskscore   : run CV-SL on placebo arm
run_cvsl_riskscore:
	Rscript code/make_checks.R 

## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects:
	Rscript code/createRDAfiles_fromSLobjects.R

## tables_figures   : create tables and figures
tables_figures:
	Rscript code/tables_figures.R

## constructSL_predict_on_vaccine   : construct SL and predict risk score on vaccinees
constructSL_predict_on_vaccine:
	Rscript code/constructSL_predict_on_vaccine.R

## get_SLweights_Modelpredictors   : get SL weights from models
get_SLweights_Modelpredictors:
	Rscript code/get_SLweights_Modelpredictors.R

## append_risk_score         : appends the computed risk scores to cleaned data
append_risk_score:
	Rscript code/append_risk_score_to_data.R

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
