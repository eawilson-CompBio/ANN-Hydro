#!/bin/bash

train_peptides=''
test_peptides=''
num_real=''
output_name=''
pred_thres=''

print_usage() {
  printf '%s\n' "Usage: bash ANN-HYDRO.sh -t peptides to use for training the network -e peptides that you wish to predict immunogencity -n number of network realizations -o names of output file -p prediction threshold [optional]"
}

while getopts 't:o:p:e:hn:' flag; do
  case "${flag}" in
    t) train_peptides="${OPTARG}" ;;
    e) test_peptides="${OPTARG}" ;;
    n) num_real="${OPTARG}" ;;
    o) output_name="${OPTARG}" ;;
    p) pred_thres="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

if [[ -z $train_peptides  ]]
then
    echo "please provided training peptides"
    exit 0
fi

if [[ -z $test_peptides  ]]
then
    echo "please provided testing peptides"
    exit 0
fi

if [[ -z $num_real  ]]
then
    echo "please specificy number of realizations"
    exit 0
fi

if [[ -z $output_name  ]]
then
    echo "please provided an output name for file"
    exit 0
fi


length_train=($(cat $train_peptides | awk -F "," 'NR>1{print $1}' | awk '{print length}' | uniq))
length_test=($(cat $test_peptides | awk -F "," 'NR>1{print $1}' | awk '{print length}' | uniq))

if [[ ${#length_train[@]} -ne 1 ]]
then
    echo "ERROR: multiple peptide lengths detected in training set. peptides need to be same length"
    exit 0
fi 

if [[ ${#length_test[@]} -ne 1 ]]
then
    echo "ERROR: multiple peptide lengths detected in testing set. peptides need to be same length"
    exit 0
fi 


if [[ $length_train -ne $length_test ]]
then
    echo "ERROR: peptides in training set are not the same length as peptides in testing set"
#    exit 0
fi

Rscript ./script/ANN_Prediction_Immunogenicity.R $train_peptides $test_peptides $num_real $output_name $pred_thres


