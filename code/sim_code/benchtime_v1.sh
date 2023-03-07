#!/bin/bash

# Activate the R Conda environment
# conda activate myRenv

# List of R scripts to run
# script_list=("npsim_tt2_tbatch.R" "nbsim_tt2_tbatch.R")

# Declare an array of sample sizes
n_list=(50 100 200 400)


# Use a for loop to print each number in the array
for number in "${n_list[@]}"
do
  Rscript npsim_compare_time_v1.R nc=1 prop_null=0.9 nsamp=$number itermax=3 filter=FALSE
done

# Loop through the script list and run each script
# for script in "${script_list[@]}"
# do

#  Rscript $script nc=4 prop_null=0.9 nsamp=50 itermax=8 filter=FALSE

#  Rscript $script nc=4 prop_null=0.9 nsamp=50 itermax=8 filter=T



# Deactivate the R Conda environment
# conda deactivate