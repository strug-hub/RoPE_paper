#!/bin/bash

# Activate the R Conda environment
# conda activate myRenv

# List of R scripts to run
script_list=("npsim_arg.R" "nbsim_arg.R")

# Loop through the script list and run each script
for script in "${script_list[@]}"
do
  echo "Running $script..."
  Rscript $script nc=8 prop_null=0.9 nsamp=50 itermax=16 fdr_control=0.05
done

# Deactivate the R Conda environment
# conda deactivate