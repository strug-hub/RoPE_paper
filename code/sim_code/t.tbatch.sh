#!/bin/bash

# Activate the R Conda environment
# conda activate myRenv

# List of R scripts to run
script_list=("npsim_tt2_tbatch.R" "nbsim_tt2_tbatch.R")

# Loop through the script list and run each script
for script in "${script_list[@]}"
do
  echo "Running FilterF $script..."
  Rscript $script nc=4 prop_null=0.9 nsamp=50 itermax=8 filter=FALSE
  echo "Running FilterT $script..."
  Rscript $script nc=4 prop_null=0.9 nsamp=50 itermax=8 filter=T

done

# Deactivate the R Conda environment
# conda deactivate