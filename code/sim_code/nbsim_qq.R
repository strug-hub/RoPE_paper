here::i_am("code/sim_code/nbsim_qq.R")

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
message(paste("This is NB simulation", sim.name))


## Simulation Settings --------------------------------------------------------
ngene <- 10000
lfc_sd <- 0.8
fdr_c.grid <- seq(0.01, 0.3, 0.01)
prop_null <- 0.8
nsamp.l <- c(50, 100, 200, 400)
filter <- T
## Load in muscle data and filter ---------------------------------------------
musc <- readRDS(here("output/tissue_data/muscle_skeletal.RDS"))
which_bad <- rowMeans(assay(musc)) < 10
musc <- musc[!which_bad, ]
epout <- readRDS(here("output/compare_powsimR/powsim_params.RDS"))
## Simulate from powsimR ------------------------------------------------------
## Get parameters from full muscle data, then simulate a smaller number of genes
## Reading the code, the first nsamp/2 individuals are in one group
## and the last nsamp/2 individuals are in the other group
res.list.NB <- list()

for (ni in 1:length(nsamp.l)) {
  nsamp <- nsamp.l[ni]
  nn <- paste0("n/2=",nsamp/2)
  set.seed(nsamp)

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
    sim.seed = nsamp,
    pLFC = function(n) rnorm(n, mean = 0, sd = lfc_sd / 2),
  )
  countdat <- psout$GeneCounts
  rownames(countdat) <- NULL
  if (any(is.na(countdat))) {
    countdat[which(is.na(countdat), arr.ind = T)] <- 0
  }

  design_mat <- cbind(1, c(rep(0, nsamp / 2), rep(1, nsamp / 2)))
  beta <- c(psout$pLFC) * 2 ## multiply back by 2 b/c divide by 2 in pLFC function

  if (filter == TRUE) {
    keep <- edgeR::filterByExpr(countdat, design = design_mat)
    countdat <- countdat[keep, ]
    beta <- beta[keep]
  }

  which_null <- abs(beta) < 10^-6

  # ## Fit methods -----------------------------------------------------
  # fitlist <- list(
  #   # vout = get_voom(countdat = countdat, design_mat = design_mat),
  #   # dout = get_DESeq2(countdat = countdat, design_mat = design_mat),
  #   # eout = get_edgeR(countdat = countdat, design_mat = design_mat),
  #   # rout = get_rope(countdat = countdat, design_mat = design_mat),
  #   # wout = get_wilcoxon(countdat = countdat, design_mat = design_mat),
  #   # nout = get_NOISeq(countdat = countdat, design_mat = design_mat)
  #   dnfout = get_DESeq2_nf(countdat = countdat, design_mat = design_mat),
  #   nfout = get_NOISeq_f(countdat = countdat, design_mat = design_mat),
  #   xout = get_dearseq(countdat = countdat, design_mat = design_mat)
  # )
  
  res <- roper::rope(datmat = countdat, X_model = design_mat, x_PI_idx = dim(design_mat)[2])
  ts <- na.omit(2 * res$adj_logLR[which_null])
  res.list.NB[[nn]] <- ts
  cat(nn)
}

saveRDS(res.list.NB, file = here("output/qq_plot/NB_qqdat.RDS"))
res.list.NP <- readRDS(here("output/qq_plot/NP_qqdat.RDS"))


par(mfrow=c(2,4))

for (ni in 1:4){
  y <- res.list.NB[[ni]]
  qqplot(qchisq(ppoints(length(y)), df = 1), y, main = (paste0("Chisq Q-Q plot for NB: n/2=",nsamp.l[ni]/2)),xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", xlim = c(0,25), ylim = c(0,25))
  ggplot2::last_plot() + qqline(y, distribution = function(p) qchisq(p, df = 1), prob = c(0.1, 0.6), col = 2)
}

for (ni in 1:4){
  y <- res.list.NP[[ni]]
  qqplot(qchisq(ppoints(length(y)), df = 1), y, main = (paste0("Chisq Q-Q plot for NP: n/2=",nsamp.l[ni]/2)),xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", xlim = c(0,25), ylim = c(0,25))
  ggplot2::last_plot() + qqline(y, distribution = function(p) qchisq(p, df = 1), prob = c(0.1, 0.6), col = 2)
}

