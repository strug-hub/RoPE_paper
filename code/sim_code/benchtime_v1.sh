#!/bin/bash

# Activate the R Conda environment
# conda activate myRenv

# Declare an array of sample sizes
n_list=(50 100 200 400)


# Use a for loop to print each number in the array
for number in "${n_list[@]}"
do
  Rscript npsim_compare_time_v1.R nc=1 prop_null=0.9 nsamp=$number itermax=3 filter=FALSE
done


# Deactivate the R Conda environment
# conda deactivate