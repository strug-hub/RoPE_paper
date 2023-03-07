here::i_am("code/figure_code/drawp03_sim_f.R")

library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggsci)
library(patchwork)

source(here("code/tools_code/draw_tool_1.R"))

## Compare fdp ----------------------------------------------------------------
# q=0.1

p.fdp.50 <- draw_fdp(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 50
)
p.fdp.100 <- draw_fdp(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 100
)
p.fdp.200 <- draw_fdp(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 200
)
p.fdp.400 <- draw_fdp(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 400
)


p1 <- p.fdp.50 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- p.fdp.100 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p.fdp.200 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p.fdp.400 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "False Discovery Proportion")
ggsave(file = here("output/sim_figures_f/fprboxp03_f.pdf"), g, width = 15, height = 5)


## Compare power ----------------------------------------------------------------
p.power.50 <- draw_power(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 50
)
p.power.100 <- draw_power(
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NB.RDS"))),
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NP.RDS"))),
  q = 0.1, n = 100
)
p.power.200 <- draw_power(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 200
)
p.power.400 <- draw_power(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NB.RDS"))),
  q = 0.1, n = 400
)


p1 <- p.power.50 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- p.power.100 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p.power.200 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p.power.400 + ylim(0, .85) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "Power")
ggsave(file = here("output/sim_figures_f/powerboxp03_f.pdf"), g, width = 15, height = 5)


# over the grid -----------------------------------------------------------

# power ---------------------------------------------------------------------
p.power.50.g <- draw_power_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NB.RDS"))),
  n = 50
)
p.power.100.g <- draw_power_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NB.RDS"))),
  n = 100
)
p.power.200.g <- draw_power_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NB.RDS"))),
  n = 200
)
p.power.400.g <- draw_power_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NB.RDS"))),
  n = 400
)


p1 <- p.power.50.g + ylim(0, .8) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- p.power.100.g + ylim(0, .8) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p.power.200.g + ylim(0, .8) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p.power.400.g + ylim(0, .8) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "Median Power", bottom = "adjusted p-value")
ggsave(file = here("output/sim_figures_f/powerp03_f.pdf"), g, width = 13, height = 8)


# FPR -------------------------------------------------------------------
p.fpr.50.g <- draw_fpr_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.50_it.48_fil.TRUE.NB.RDS"))),
  n = 50
)
p.fpr.100.g <- draw_fpr_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.100_it.48_fil.TRUE.NB.RDS"))),
  n = 100
)
p.fpr.200.g <- draw_fpr_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.200_it.48_fil.TRUE.NB.RDS"))),
  n = 200
)
p.fpr.400.g <- draw_fpr_grid(
  NP.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NP.RDS"))),
  NB.df = as_tibble(readRDS(file = here("output/sim_df_out_f/nc.8_pn0.7_ns.400_it.48_fil.TRUE.NB.RDS"))),
  n = 400
)


p1 <- p.fpr.50.g + ylim(0, .7) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p2 <- p.fpr.100.g + ylim(0, .7) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p.fpr.200.g + ylim(0, .7) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p.fpr.400.g + ylim(0, .7) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "Median false discovery proportion", bottom = "adjusted p-value")
ggsave(file = here("output/sim_figures_f/fprp03_f.pdf"), g, width = 13, height = 8)
