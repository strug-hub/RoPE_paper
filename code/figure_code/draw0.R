here::i_am("code/figure_code/draw0.R")

library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
powsimr_df <- as_tibble(readRDS(file = here("output/diff_exp_out/nc.4_pn0.9_ns.50_it.8_q0.05.NB.RDS")))
seqgendiff_df <- as_tibble(readRDS(file = here("output/diff_exp_out/nc.4_pn0.9_ns.50_it.8_q0.05.NP.RDS")))

powsimr_df %>%
  group_by(seed) %>%
  count() %>%
  filter(n > 1) %>%
  nrow() ->
  numnonunique
stopifnot(numnonunique == 0)

## Compare median PVE on non-null genes ---------------------------------------
powsimr_df %>%
  select(lfc_sd, mpve) %>%
  mutate(method = "powsimR") ->
  pstemp


seqgendiff_df %>%
  select(mpve) %>%
  mutate(method = "seqgendiff", lfc_sd = 0.8) %>%
  select(lfc_sd, mpve, method) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(method_lfc = method) %>%
  ggplot(aes(x = method_lfc, y = mpve)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Simulation Method") +
  ylab("Median PVE") ->
  pl

pl

## Compare fdp ----------------------------------------------------------------
powsimr_df %>%
  select(contains("fpr"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "fdp") %>%
  mutate(method = str_replace(method, "fpr_", ""),
         sim = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(contains("fpr")) %>%
  gather(key = "method", value = "fdp") %>%
  mutate(method = str_replace(method, "fpr_", ""),
         sim = "seqgendiff",
         lfc_sd = 0.8) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = sim) %>%
  mutate(method = recode(method,
                         "dout" = "DESeq2",
                         "eout" = "edgeR",
                         "vout" = "voom+limma",
                         "rout" = "RoPE",
                         "wout" = "Wilcoxon",
                         "nout" = "NOISeq")) %>%
  ggplot(aes(x = method, y = fdp, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  geom_hline(yintercept = 0.05, lty = 2) +
  xlab("Method") +
  ylab("FDP") ->
  pl

pl

## Compare power ----------------------------------------------------------------
powsimr_df %>%
  select(contains("power"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "power") %>%
  mutate(method = str_replace(method, "power_", ""),
         sim = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(contains("power")) %>%
  gather(key = "method", value = "power") %>%
  mutate(method = str_replace(method, "power_", ""),
         sim = "seqgendiff",
         lfc_sd = 0.8) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = sim) %>%
  mutate(method = recode(method,
                         "dout" = "DESeq2",
                         "eout" = "edgeR",
                         "vout" = "voom+limma",
                         "rout" = "RoPE",
                         "wout" = "Wilcoxon",
                         "nout" = "NOISeq")) %>%
  ggplot(aes(x = method, y = power, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  xlab("Method") +
  ylab("Power") ->
  pl

pl
