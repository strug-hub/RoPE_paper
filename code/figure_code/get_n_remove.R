here::i_am("code/figure_code/get_n_remove.R")

library(here)
suppressPackageStartupMessages(library(tidyverse))

source(here("code/tools_code/n.remove.R"))

# non-filtered ----------------------------------------------------------------


n.r.50 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.50_it.48_fil.FALSE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.50_it.48_fil.FALSE.NB.RDS"))),
  n = 50
)
n.r.100 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.100_it.48_fil.FALSE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.100_it.48_fil.FALSE.NB.RDS"))),
  n = 100
)
n.r.200 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.200_it.48_fil.FALSE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.200_it.48_fil.FALSE.NB.RDS"))),
  n = 200
)
n.r.400 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.400_it.48_fil.FALSE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out/nc.8_pn0.7_ns.400_it.48_fil.FALSE.NB.RDS"))),
  n = 400
)

# filtered ----------------------------------------------------------------


n.r.f50 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NB.RDS"))),
  n = 50
)
n.r.f100 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NB.RDS"))),
  n = 100
)
n.r.f200 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NB.RDS"))),
  n = 200
)
n.r.f400 <- get.n.remove(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NB.RDS"))),
  n = 400
)

bind_rows(
  n.r.50, 
  n.r.100, 
  n.r.200,
  n.r.400,
) -> N.remove.nf

bind_rows(
  n.r.f50, 
  n.r.f100, 
  n.r.f200,
  n.r.f400,
) -> N.remove.f

N.remove.nf %>% pivot_wider(names_from = size, values_from = m_remove)
N.remove.f %>% pivot_wider(names_from = size, values_from = m_remove)
