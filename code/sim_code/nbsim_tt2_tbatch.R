here::i_am("code/sim_code/nbsim_tt2.R")

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
message(paste("This is NB simulation",sim.name))


## Simulation Settings --------------------------------------------------------
ngene <- 10000
lfc_sd <- 0.8
fdr_c.grid <- seq(0.01, 0.3, 0.01)
pardf <- expand.grid(seed = seq_len(itermax), lfc_sd = lfc_sd)

## Load in muscle data and filter ---------------------------------------------
musc <- readRDS(here("output/tissue_data/muscle_skeletal.RDS"))
which_bad <- rowMeans(assay(musc)) < 10
musc <- musc[!which_bad, ]

## Simulate from powsimR ------------------------------------------------------
## Get parameters from full muscle data, then simulate a smaller number of genes
## Reading the code, the first nsamp/2 individuals are in one group
## and the last nsamp/2 individuals are in the other group
epout <- readRDS(here("output/compare_powsimR/powsim_params.RDS"))

cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

retmat <- foreach(
  iterindex = seq_len(nrow(pardf)),
  .combine = rbind,
  .export = c("simulateCounts")
) %dopar% {
  set.seed(pardf[iterindex, "seed"])
  lfc_sd <- pardf[iterindex, "lfc_sd"]

  ## Simulate data ---------------------------------------------------
  ## We divide by lfc_sd by 2 in powsimR but not in seqgendiff because
  ## the design matrix (used to simulate counts) in powsimR is c(rep(-1, nsamp/2), rep(1, nsamp/2))
  ## while the design matrix in seqgendiff if c(rep(0, nsamp/2), rep(1, nsamp/2)).
  ## Dividing by 2 will allow us to use the 0/1 design matrix in voom-limma
  ## rather than the -1/1 design matrix.
  psout <- simulateCounts(
    n = c(nsamp / 2, nsamp / 2),
    ngenes = ngene,
    p.DE = 1 - prop_null,
    params = epout,
    sim.seed = iterindex,
    pLFC = function(n) rnorm(n, mean = 0, sd = lfc_sd / 2),
  )
  countdat <- psout$GeneCounts
  rownames(countdat) <- NULL
  if (any(is.na(countdat))) {
    countdat[which(is.na(countdat), arr.ind = T)] <- 0
  }

  design_mat <- cbind(1, c(rep(0, nsamp / 2), rep(1, nsamp / 2)))
  beta <- c(psout$pLFC) * 2 ## multiply back by 2 b/c divide by 2 in pLFC function
  
  if (filter == TRUE){
    keep <- edgeR::filterByExpr(countdat, design = design_mat)
    countdat <- countdat[keep,]
    beta <- beta[keep]
  }
  
  which_null <- abs(beta) < 10^-6

  # ## Fit methods -----------------------------------------------------
  # fitlist <- list(
  #   vout = get_voom(countdat = countdat, design_mat = design_mat),
  #   dout = get_DESeq2(countdat = countdat, design_mat = design_mat),
  #   eout = get_edgeR(countdat = countdat, design_mat = design_mat),
  #   rout = get_rope(countdat = countdat, design_mat = design_mat),
  #   wout = get_wilcoxon(countdat = countdat, design_mat = design_mat)
  #   # nout = get_NOISeq(countdat = countdat, design_mat = design_mat)
  # )
  # 
  # ## Assess fits -----------------------------------------------------
  # res.perf <- unlist(lapply(fdr_c.grid, FUN = function(fdr_c) {
  #   fitlist <- lapply(fitlist, FUN = function(obj) {
  #     # obj$qval <- p.adjust(obj$pval, method = "BH")
  #     obj$discovery <- obj$qval < fdr_c
  #     return(obj)
  #   })
  # 
  #   fprvec <- sapply(fitlist, FUN = function(obj) {
  #     if (any(obj$discovery, na.rm = TRUE)) {
  #       mean(which_null[obj$discovery], na.rm = TRUE)
  #     } else {
  #       0
  #     }
  #   })
  #   names(fprvec) <- paste0("fpr_", names(fprvec), "_q_", fdr_c)
  # 
  #   powervec <- sapply(fitlist, FUN = function(obj) {
  #     mean(obj$discovery[!which_null], na.rm = TRUE)
  #   })
  #   names(powervec) <- paste0("power_", names(powervec), "_q_", fdr_c)
  #   return(c(fprvec, powervec))
  # }))
  # 
  # msevec <- sapply(fitlist, FUN = function(obj) {
  #   mean((obj$bhat - beta)^2, na.rm = TRUE)
  # })
  # names(msevec) <- paste0("mse_", names(msevec))
  # 
  # ## Summary stat of count matrix ------------------------------------
  # varvec <- apply(log2(countdat + 0.5)[!which_null, , drop = FALSE], 1, var)
  # betavarvec <- apply(tcrossprod(beta[!which_null], design_mat[, 2]), 1, var)
  # mpve <- median(betavarvec / varvec)
  # 
  # ## Return ----------------------------------------------------------
  # retvec <- c(res.perf ,msevec, mpve = mpve, unlist(pardf[iterindex, ]))

  retvec <- MatrixGenerics::colSums2(countdat)
  
  retvec
}
stopCluster(cl)

saveRDS(object = retmat, file = here("output", "batch_test", paste0(sim.name, "fake.NB.RDS")))
