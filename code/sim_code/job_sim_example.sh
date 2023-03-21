#!/bin/bash
#SBATCH -N 1 -c 9 
#SBATCH --mem 100G
#SBATCH --time 48:00:00 
#SBATCH --job-name a_sim_f8

# Activate the R Conda environment
eval "$(conda shell.bash hook)"
conda activate myRenv

# List of R scripts to run
script_list=("npsim_code.R" "nbsim_code.R")

# Declare an array of sample sizes
n_list=(50 100 200 400)

# Use nested for loops to print all combinations
for number in "${n_list[@]}"
do
  for script in "${script_list[@]}"
  do
    echo "Running $script $number..."
	# variables: nc: num of cores; prop_null: proportion of non-DE genes; itermax: number of replicates per setting; filter: low-count filter
    Rscript $script nc=8 prop_null=0.8 nsamp=$number itermax=48 filter=T
  done
done

echo "Job complete, captain!" 

# Deactivate the R Conda environment
conda deactivate