here::i_am("code/sim_code/npsim_qq.R")

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
prop_null <- 0.8
nsamp.l <- c(50, 100, 200, 400)
filter <- T
## Load in muscle data and filter ---------------------------------------------
musc <- readRDS(here("output/tissue_data/muscle_skeletal.RDS"))
which_bad <- rowMeans(assay(musc)) < 10
musc <- musc[!which_bad, ]
fullcounts <- assay(musc)

res.list <- list()

for (ni in 1:length(nsamp.l)) {
  nsamp <- nsamp.l[ni]
  nn <- paste0("n/2=",nsamp/2)
  set.seed(nsamp)
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
  res <- roper::rope(datmat = countdat, X_model = design_mat, x_PI_idx = dim(design_mat)[2])
  ts <- na.omit(2 * res$adj_logLR[which_null])
  res.list[[nn]] <- ts
  cat(nn)
}

saveRDS(res.list, file = here("output/qq_plot/NP_qqdat.RDS"))

res.list$`n/2=25`
y <- res.list$`n/2=200`

par(mfrow=c(1,4))
for (ni in 1:4){
  y <- res.list[[ni]]
  qqplot(qchisq(ppoints(length(y)), df = 1), y, main = (paste0("Chisq Q-Q plot for NP: n/2=",nsamp.l[ni]/2)),xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", xlim = c(0,25), ylim = c(0,25))
  ggplot2::last_plot() + qqline(y, distribution = function(p) qchisq(p, df = 1), prob = c(0.1, 0.6), col = 2)
}

qqplot(qchisq(ppoints(length(y)), df = 1), y, main = expression("Q-Q plot for" ~~ {chi^2}[nu == 1]), )

library(ggplot2)
dev.off()
