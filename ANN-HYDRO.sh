#!/bin/bash



train_peptides=''
test_peptides=''
num_real=''
output_name=''
pred_thres=''
gui=''

print_usage() {
  printf '%s\n' "Usage: bash ANN-HYDRO.sh -t peptides to use for training the network -e peptides that you wish to predict immunogencity -n number of network realizations -o names of output file -p prediction threshold [optional] -g gui interface"
}

while getopts 't:o:p:e:hn:g:' flag; do
  case "${flag}" in
    t) train_peptides="${OPTARG}" ;;
    e) test_peptides="${OPTARG}" ;;
    n) num_real="${OPTARG}" ;;
    o) output_name="${OPTARG}" ;;
    p) pred_thres="${OPTARG}" ;;
    g) gui="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done


if [[ -z $gui ]]
then
    

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

 if [[ -z $pred_thres  ]]
  then
      echo "no prediction threshold provided. setting to 0.5"
      pred_thres=0.5
 fi

 if [[ $pred_thres -gt 1 ]]
 then
      echo "prediction threshold must be a number between 0 and 1. setting to 0.5"
      pred_thres=0.5
 elif [[ $pred_thres -lt 0 ]]
 then
     
      echo "prediction threshold must be a number between 0 and 1. setting to 0.5"
      pred_thres=0.5
 fi
      
      
  
else
    zenity --info --text="select peptide training set" --width 150 --height 150 
	   train_peptides=$(zenity --file-selection --filename=$PWD)
    zenity --info --text="select peptide testing set" --width 150 --height 150
    test_peptides=$(zenity --file-selection --title="select peptide testing set"--filename=$PWD)
    param=($(zenity --forms --title="prediction parameters" \
	   --separator=" " \
	   --add-entry="number of realizations (e.g 30)" \
	   --add-entry="output file name (e.g HLA-A2_predicted.txt)" \
	   --add-entry="prediction threshold (e.g 0.5 )" \
	   ))
       if [[ -z ${param[0]}  ]]
       then
	    
	   echo "please specificy number of realizations"
	   exit 0
  
       fi

       if [[ -z ${param[1]}  ]]
       then
	   echo "please provided an output name for file"
	   exit 0
       fi

       if [[ -z ${param[2]}  ]]
       then
	   echo "no prediction threshold provided. setting to 0.5"
	   pred_thres=0.5
       fi

       if [[ $pred_thres -gt 1 ]]
       then
	   echo "prediction threshold must be a number between 0 and 1. setting to 0.5"
	   pred_thres=0.5
       elif [[ $pred_thres -lt 0 ]]
       then
	   
	    echo "prediction threshold must be a number between 0 and 1. setting to 0.5"
	    pred_thres=0.5
       fi
   
    num_real=${param[0]}
    output_name=${param[1]}
    pred_thres=${param[2]}
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

dir=$(greadlink -f ./Predictions/)
zenity --info --text=$(echo "prediction complete and saved to"$dir)
