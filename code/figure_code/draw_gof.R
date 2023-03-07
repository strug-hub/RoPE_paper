here::i_am("code/figure_code/draw_gof.R")

library(here)
library(tidyverse)

here("output/gof_f/nc.4_pn0.8_ns.50_it.48_fil.TRUE.NB.RDS")
here("output/gof_f/nc.4_pn0.8_ns.50_it.48_fil.TRUE.NP.RDS")


NB.list <- readRDS(here("output/gof_f/nc.4_pn0.8_ns.50_it.48_fil.TRUE.NB.RDS"))
NP.list <- readRDS(here("output/gof_f/nc.4_pn0.8_ns.50_it.48_fil.TRUE.NP.RDS"))

obj <- NB.list[[1]]



-log10(obj$null.p$g0)

-log10(obj$null.p$g1)

