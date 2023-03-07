here::i_am("code/sim_code/get_gtex.R")

library(here)
library(tidyverse)

# recount3 ----------------------------------------------------------------
library(recount3)

human_projects <- available_projects()

# dim(human_projects)
# head(human_projects)
#
# subset(human_projects, file_source == "gtex" & project_type == "data_sources")

proj_info <- subset(
  human_projects,
  project == "MUSCLE" & project_type == "data_sources"
)

rse_gene_gtex_muscle<- create_rse(proj_info)
rse_gene_gtex_muscle

assay(rse_gene_gtex_muscle, "counts") <- transform_counts(rse_gene_gtex_muscle)


# create base count -------------------------------------------------------
# library(edgeR)

# musc <- assay(rse_gene_gtex_muscle, "counts")
musc <- rse_gene_gtex_muscle
# which_bad <- rowMeans(assay(rse_gene_gtex_muscle)) < 10
# dat.count <- dat.count.0[!which_bad, ]
# rm(dat.count.0)
# or
# keep <- filterByExpr(dat.count.0)
# musc <- dat.count.0[keep, ]
# rm(dat.count.0)
#
saveRDS(object = musc, file = here("output/tissue_data/muscle_skeletal.RDS"))