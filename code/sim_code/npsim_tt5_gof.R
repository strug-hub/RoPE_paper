here::i_am("code/sim_code/npsim_tt5_gof.R")

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
} else if (length(args) == 5) {
  eval(parse(text = args[[1]]))
  eval(parse(text = args[[2]]))
  eval(parse(text = args[[3]]))
  eval(parse(text = args[[4]]))
  eval(parse(text = args[[5]]))
} else {
  stop("This code needs to be provided with either 0 or 5 arg", call. = F)
}

sim.name <- paste0("nc.", nc, "_pn", prop_null, "_ns.", nsamp, "_it.", itermax, "_fil.", filter)
message(paste("This is NP simulation", sim.name))

## Simulation Settings --------------------------------------------------------
# prop_null <- 0.9
# nsamp     <- 50
# itermax   <- 500
# itermax   <- 8
# fdr_control <- 0.05
ngene <- 10000
lfc_sd <- 0.8
fdr_c.grid <- seq(0.01, 0.3, 0.01)
## Load in muscle data and filter ---------------------------------------------
musc <- readRDS(here("output/tissue_data/muscle_skeletal.RDS"))
which_bad <- rowMeans(assay(musc)) < 10
musc <- musc[!which_bad, ]
fullcounts <- assay(musc)

## Set up parallel computing environment --------------------------------------
# nc <- 4
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

retmat <- foreach(
  iterindex = seq_len(itermax)
) %dopar% {
  set.seed(iterindex)
  ## Simulate data ---------------------------------------------------
  which_gene <- sort(sample(seq_len(BiocGenerics::nrow(musc)), ngene))
  which_samp <- sort(sample(seq_len(BiocGenerics::ncol(musc)), nsamp))
  submusc <- fullcounts[which_gene, which_samp]
  thout <- seqgendiff::thin_2group(
    mat = submusc,
    prop_null = prop_null,
    signal_fun = stats::rnorm,
    signal_params = list(mean = 0, sd = lfc_sd),
    group_prop = 0.5
  )
  countdat <- thout$mat
  design_mat <- cbind(thout$design_obs, thout$designmat)
  beta <- c(thout$coefmat)

  if (filter == TRUE) {
    keep <- edgeR::filterByExpr(countdat, design = design_mat)
    countdat <- countdat[keep, ]
    beta <- beta[keep]
  }

  which_null <- abs(beta) < 10^-6
  dge <- edgeR::DGEList(counts = countdat, group = design_mat[, 2])
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(y = dge, design = design_mat)
  
  count.norm.0 <- edgeR::cpm(dge)[,design_mat[, 2] == 0]
  count.norm.1 <- edgeR::cpm(dge)[,design_mat[, 2] == 1]
  
  dis <- dge$common.dispersion
  # dis <- dge$tagwise.dispersion
  
  library(vcd)
  
  p.0 <- unlist(lapply(1:nrow(countdat), FUN = function(gene) {
    tt <- summary(goodfit(round(count.norm.0)[gene,], type = "nbinomial", par = list(size = 1/dis)))
    tt[,3]}))
  
  p.1 <- unlist(lapply(1:nrow(countdat), FUN = function(gene) {
    tt <- summary(goodfit(round(count.norm.1)[gene,], type = "nbinomial", par = list(size = 1/dis)))
    tt[,3]}))
  
  ## Return ----------------------------------------------------------
  null.p <- tibble(g0 = p.0[which_null],g1 = p.1[which_null])
  nnull.p <- tibble(g0 = p.0[!which_null],g1 = p.1[!which_null])
  res <- list(null.p = null.p, nnull.p= nnull.p)
  res
}
stopCluster(cl)

s.folder.name <- "gof"

if (filter == FALSE) {
  s.path <- here("output", s.folder.name, paste0(sim.name, ".NP.RDS"))
} else {
  s.path <- here("output", paste0(s.folder.name, "_f"), paste0(sim.name, ".NP.RDS"))
}

saveRDS(object = retmat, file = s.path)
