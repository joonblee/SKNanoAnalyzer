#!/bin/bash

##############
### setups ###
##############

#Eras=(2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix)
Eras=(2022 2022EE 2023 2023BPix)
Eras=(2022EE 2023 2023BPix)

Sample_Run2="SampleLists/Run2mc_all.txt"
Sample_Run3="SampleLists/Run3mc.txt"
Sample_Run3="SampleLists/Run3signal.txt"

nBatch=-1

list_to_input() {
  local listfile="$1"
  awk '
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      print $0
    }
  ' "$listfile" | paste -sd, -
}

for era in "${Eras[@]}"
do
  if [[ "$era" == "2016preVFP" || "$era" == "2016postVFP" || "$era" == "2017" || "$era" == "2018" ]]; then
    Sample="$Sample_Run2"
    DataSample="SingleMuon"
  elif [[ "$era" == "2022" || "$era" == "2022EE" ]]; then
    Sample="$Sample_Run3"
    DataSample="Muon"
  elif [[ "$era" == "2023" || "$era" == "2023BPix" ]]; then
    Sample="$Sample_Run3"
    DataSample="Muon0,Muon1"
  fi

  mc_samples=$(list_to_input "$Sample")
  #input_samples="${DataSample},${mc_samples}"
  input_samples="${mc_samples}"

  SKNano.py \
    -a Skim_NIsoMuon \
    -i "$input_samples" \
    -n "$nBatch" \
    -e "$era" \
    --nmax 1000 \
    --skimming_mode #\
    #&> "submit_skim_${era}.log" 

  sleep 2
done
