here::i_am("code/figure_code/draw1.R")

library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggsci)
powsimr_df <- as_tibble(readRDS(file = here("output/diff_exp_out/nc.8_pn0.9_ns.50_it.48_q0.05.NB.RDS")))
seqgendiff_df <- as_tibble(readRDS(file = here("output/diff_exp_out/nc.8_pn0.9_ns.50_it.48_q0.05.NP.RDS")))

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
# q=0.1
powsimr_df %>%
  select(matches("^fpr.*_q_0\\.1$"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "fdp") %>%
  mutate(
    method = str_replace(method, "fpr_", ""),
    sim = "powsimR"
  ) ->
pstemp

seqgendiff_df %>%
  select(matches("^fpr.*_q_0\\.1$")) %>%
  gather(key = "method", value = "fdp") %>%
  mutate(
    method = str_replace(method, "fpr_", ""),
    sim = "seqgendiff",
    lfc_sd = 0.8
  ) ->
sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = sim) %>%
  mutate(method = recode(method,
    "dout_q_0.1" = "DESeq2",
    "eout_q_0.1" = "edgeR",
    "vout_q_0.1" = "voom+limma",
    "rout_q_0.1" = "RoPE",
    "wout_q_0.1" = "Wilcoxon",
    "nout_q_0.1" = "NOISeq"
  )) %>%
  ggplot(aes(x = method, y = fdp, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  geom_hline(yintercept = 0.1, lty = 2) +
  xlab("Method") +
  ylab("FDP") ->
pl

pl

## Compare power ----------------------------------------------------------------
powsimr_df %>%
  select(matches("^power_.*_q_0\\.1$"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "power") %>%
  mutate(
    method = str_replace(method, "power_", ""),
    sim = "powsimR"
  ) ->
pstemp

seqgendiff_df %>%
  select(matches("^power_.*_q_0\\.1$")) %>%
  gather(key = "method", value = "power") %>%
  mutate(
    method = str_replace(method, "power_", ""),
    sim = "seqgendiff",
    lfc_sd = 0.8
  ) ->
sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = sim) %>%
  mutate(method = recode(method,
    "dout_q_0.1" = "DESeq2",
    "eout_q_0.1" = "edgeR",
    "vout_q_0.1" = "voom+limma",
    "rout_q_0.1" = "RoPE",
    "wout_q_0.1" = "Wilcoxon",
    "nout_q_0.1" = "NOISeq"
  )) %>%
  ggplot(aes(x = method, y = power, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  xlab("Method") +
  ylab("Power") ->
pl

pl


# over the grid -----------------------------------------------------------

# power ---------------------------------------------------------------------
tt.nb.power <- powsimr_df %>%
  select(matches("^power_.*")) %>% mutate(iter = 1:nrow(powsimr_df)) %>%
  select(matches("^power_.*"),iter) %>%
  pivot_longer(cols = !iter, names_to = "column_name", values_to = "power") %>%
  separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
  mutate(method = gsub("power_", "", method)) %>%
  mutate(method = factor(recode(method,
                                "dout" = "DESeq2",
                                "eout" = "edgeR",
                                "vout" = "voom+limma",
                                "rout" = "RoPE",
                                "wout" = "Wilcoxon",
                                "nout" = "NOISeq"
  ))) %>%
  mutate(method = fct_relevel(
    method,
    "RoPE",
    "edgeR",
    "DESeq2",
    "voom+limma",
    "NOISeq",
    "Wilcoxon"
  ))

power.nb <- tt.nb.power %>% group_by(method, qvalue) %>% summarise(med_power = median(power)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue))


tt.np.power <- seqgendiff_df %>%
  select(matches("^power_.*")) %>% mutate(iter = 1:nrow(powsimr_df)) %>%
  select(matches("^power_.*"),iter) %>%
  pivot_longer(cols = !iter, names_to = "column_name", values_to = "power") %>%
  separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
  mutate(method = gsub("power_", "", method)) %>%
  mutate(method = factor(recode(method,
                                "dout" = "DESeq2",
                                "eout" = "edgeR",
                                "vout" = "voom+limma",
                                "rout" = "RoPE",
                                "wout" = "Wilcoxon",
                                "nout" = "NOISeq"
  ))) %>%
  mutate(method = fct_relevel(
    method,
    "RoPE",
    "edgeR",
    "DESeq2",
    "voom+limma",
    "NOISeq",
    "Wilcoxon"
  ))

power.np <- tt.np.power %>% group_by(method, qvalue) %>% summarise(med_power = median(power)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue))

power.dat <- bind_rows(power.nb %>% mutate(sim = "NB-parametric"),
                       power.np %>% mutate(sim = "Non-parametric"))

ggplot(power.dat, aes(x = qvalue, y = med_power, color = method)) + geom_line(linewidth = 0.8) + theme_bw()  + 
  facet_grid(sim~.) + scale_color_npg()


# FDR -------------------------------------------------------------------

tt.nb <- powsimr_df %>%
  select(matches("^fpr_.*")) %>% mutate(iter = 1:nrow(powsimr_df)) %>%
  select(matches("^fpr_.*"),iter) %>%
  pivot_longer(cols = !iter, names_to = "column_name", values_to = "fpr") %>%
  separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
  mutate(method = gsub("fpr_", "", method)) %>%
  mutate(method = factor(recode(method,
                                "dout" = "DESeq2",
                                "eout" = "edgeR",
                                "vout" = "voom+limma",
                                "rout" = "RoPE",
                                "wout" = "Wilcoxon",
                                "nout" = "NOISeq"
  ))) %>%
  mutate(method = fct_relevel(
    method,
    "RoPE",
    "edgeR",
    "DESeq2",
    "voom+limma",
    "NOISeq",
    "Wilcoxon"
  ))

fpr.nb <- tt.nb %>% group_by(method, qvalue) %>% summarise(med_fpr = median(fpr)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue))


tt.np <- seqgendiff_df %>%
  select(matches("^fpr_.*")) %>% mutate(iter = 1:nrow(powsimr_df)) %>%
  select(matches("^fpr_.*"),iter) %>%
  pivot_longer(cols = !iter, names_to = "column_name", values_to = "fpr") %>%
  separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
  mutate(method = gsub("fpr_", "", method)) %>%
  mutate(method = factor(recode(method,
                                "dout" = "DESeq2",
                                "eout" = "edgeR",
                                "vout" = "voom+limma",
                                "rout" = "RoPE",
                                "wout" = "Wilcoxon",
                                "nout" = "NOISeq"
  ))) %>%
  mutate(method = fct_relevel(
    method,
    "RoPE",
    "edgeR",
    "DESeq2",
    "voom+limma",
    "NOISeq",
    "Wilcoxon"
  ))

fpr.np <- tt.np %>% group_by(method, qvalue) %>% summarise(med_fpr = median(fpr)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue))

fpr.dat <- bind_rows(fpr.nb %>% mutate(sim = "NB-parametric"),
                     fpr.np %>% mutate(sim = "Non-parametric"))

ggplot(fpr.dat, aes(x = qvalue, y = med_fpr, color = method)) + geom_line(linewidth = 0.8) + theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  facet_grid(sim~.) + scale_color_npg()

