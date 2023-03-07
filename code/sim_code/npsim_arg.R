here::i_am("code/sim_code/npsim_arg.R")

library(here)
suppressPackageStartupMessages(library(tidyverse))

##########################
## Simulate data in two-group model, compare basic approaches
##########################

## load packages --------------------------------------------------------------
library(seqgendiff)
suppressPackageStartupMessages(library(SummarizedExperiment))
source(here("code/tools_code/de_methods_rope.R"))
suppressPackageStartupMessages(library(doSNOW))

# Number of threads to use for multithreaded computing. This must be
# specified in the command-line shell; e.g., to use 8 threads, run
# command
#
#  R CMD BATCH '--args nc=8' mouthwash_sims.R
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
  prop_null <- 0.9
  nsamp <- 50
  itermax <- 8
  # fdr_control <- 0.05
  filter <- FALSE
} else if (length(args) == 5){
  eval(parse(text = args[[1]]))
  eval(parse(text = args[[2]]))
  eval(parse(text = args[[3]]))
  eval(parse(text = args[[4]]))
  eval(parse(text = args[[5]]))
} else {
  stop("This code needs to be provided with either 0 or 5 arg",call. = F)
}

sim.name <- paste0("nc.", nc, "_pn", prop_null, "_ns.", nsamp, "_it.", itermax,"_fil.",filter)
message(paste("This is NP simulation",sim.name))
