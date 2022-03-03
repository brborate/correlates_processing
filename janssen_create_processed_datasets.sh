#!/bin/bash
declare -a arr=("janssen_pooled_real" "janssen_pooled_realADCP" 
                "janssen_na_real" "janssen_na_realADCP" 
                "janssen_la_real" "janssen_la_realADCP" 
                "janssen_sa_real" "janssen_sa_realADCP")

for i in "${arr[@]}"
do
   echo "$i"
   export TRIAL=$i
   make data_processed
   make deploy_processed_dataset
   echo "Analysis dataset created for $i and deployed!"
done