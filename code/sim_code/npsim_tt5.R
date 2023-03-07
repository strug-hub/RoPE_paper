here::i_am("code/sim_code/npsim_tt5.R")

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
  iterindex = seq_len(itermax),
  .combine = rbind
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

  ## Fit methods -----------------------------------------------------
  fitlist <- list(
    # vout = get_voom(countdat = countdat, design_mat = design_mat),
    # dout = get_DESeq2(countdat = countdat, design_mat = design_mat),
    # eout = get_edgeR(countdat = countdat, design_mat = design_mat),
    # rout = get_rope(countdat = countdat, design_mat = design_mat),
    # wout = get_wilcoxon(countdat = countdat, design_mat = design_mat),
    # nout = get_NOISeq(countdat = countdat, design_mat = design_mat)
    dnfout = get_DESeq2_nf(countdat = countdat, design_mat = design_mat),
    nfout = get_NOISeq_f(countdat = countdat, design_mat = design_mat),
    xout = get_dearseq(countdat = countdat, design_mat = design_mat)
  )

  ## Assess fits -----------------------------------------------------
  res.perf <- unlist(lapply(fdr_c.grid, FUN = function(fdr_c) {
    fitlist <- lapply(fitlist, FUN = function(obj) {
      # obj$qval <- p.adjust(obj$pval, method = "BH")
      obj$discovery <- obj$qval < fdr_c
      return(obj)
    })

    fprvec <- sapply(fitlist, FUN = function(obj) {
      if (any(obj$discovery, na.rm = TRUE)) {
        mean(which_null[obj$discovery], na.rm = TRUE)
      } else {
        0
      }
    })
    names(fprvec) <- paste0("fpr_", names(fprvec), "_q_", fdr_c)

    powervec <- sapply(fitlist, FUN = function(obj) {
      mean(obj$discovery[!which_null], na.rm = TRUE)
    })
    names(powervec) <- paste0("power_", names(powervec), "_q_", fdr_c)
    return(c(fprvec, powervec))
  }))


  msevec <- sapply(fitlist, FUN = function(obj) {
    suppressWarnings(mean((obj$bhat - beta)^2, na.rm = TRUE))
  })
  names(msevec) <- paste0("mse_", names(msevec))

  # Assess outliers removal:
  olnavec <- sapply(fitlist, FUN = function(obj) {
    as.integer(obj$n.remove)
  })
  names(olnavec) <- paste0("olna_", names(olnavec))

  ## Summary stat of count matrix ------------------------------------
  varvec <- apply(log2(countdat + 0.5)[!which_null, , drop = FALSE], 1, var)
  betavarvec <- apply(tcrossprod(beta[!which_null], design_mat[, 2]), 1, var)
  mpve <- median(betavarvec / varvec)


  ## Return ----------------------------------------------------------
  retvec <- c(res.perf, msevec, mpve = mpve, olnavec)

  retvec
}
stopCluster(cl)

s.folder.name <- "sim_ad_out"

if (filter == FALSE) {
  s.path <- here("output", s.folder.name, paste0(sim.name, ".NP.RDS"))
} else {
  s.path <- here("output", paste0(s.folder.name, "_f"), paste0(sim.name, ".NP.RDS"))
}

saveRDS(object = retmat, file = s.path)
