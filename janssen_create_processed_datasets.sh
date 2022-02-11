#!/bin/bash
declare -a arr=("janssen_pooled_realbAb" "janssen_pooled_realADCP" "janssen_pooled_realPsV"
                "janssen_na_realbAb" "janssen_na_realADCP" "janssen_na_realPsV"
                "janssen_la_realbAb" "janssen_la_realADCP" "janssen_la_realPsV"
                "janssen_sa_realbAb" "janssen_sa_realADCP" "janssen_sa_realPsV")

for i in "${arr[@]}"
do
   echo "$i"
   export TRIAL=$i
   make data_processed
   make deploy_processed_dataset
   echo "Analysis dataset created for $i and deployed!"
done