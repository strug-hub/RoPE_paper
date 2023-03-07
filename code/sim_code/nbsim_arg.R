here::i_am("code/sim_code/nbsim_tt0.R")

library(here)
suppressPackageStartupMessages(library(tidyverse))
source(here("code/tools_code/Simulate.R"))
source(here("code/tools_code/utils_misc.R"))
source(here("code/tools_code/utils_sim.R"))
source(here("code/tools_code/de_methods_rope.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(doSNOW))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
  prop_null <- 0.9
  nsamp <- 50
  itermax <- 8
  fdr_control <- 0.05
} else {
  eval(parse(text = args[[1]]))
  eval(parse(text = args[[2]]))
  eval(parse(text = args[[3]]))
  eval(parse(text = args[[4]]))
  eval(parse(text = args[[5]]))
}

sim.name <- paste0("nc.", nc, "_pn", prop_null, "_ns.", nsamp, "_it.", itermax, "_q", fdr_control)
message(sim.name)
