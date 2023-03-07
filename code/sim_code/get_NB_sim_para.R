here::i_am("code/sim_code/get_NB_sim_para.R")
library(here)

# ## Estimate parameters for powsimR simulations --------------------------
# suppressPackageStartupMessages(library(SummarizedExperiment))
musc <- readRDS(here("output/tissue_data/muscle_skeletal.RDS"))
source(here('code/tools_code/Simulate.R'))
source(here('code/tools_code/utils_misc.R'))
source(here('code/tools_code/utils_sim.R'))
source(here('code/tools_code/utils_param.R'))
source(here('code/tools_code/utils_normalise.R'))
source(here('code/tools_code/Param.R'))

which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
epout <- estimateParam(countData = assay(musc),
                       Distribution = "NB",
                       normalisation = "TMM",
                       RNAseq = "bulk")

message("NB parameters estimation is complete.")

saveRDS(epout, file = here("output/compare_powsimR/powsim_params.RDS"))
