library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggsci)

get.n.remove <- function(NP.df, NB.df, n = 50) {
  NB.df %>%
    select(matches("^olna_.*"), lfc_sd) %>%
    gather(-lfc_sd, key = "method", value = "n_remove") %>%
    mutate(
      method = str_replace(method, "olna_", ""),
      sim = "NB-parametric"
    ) ->
  pstemp


  NP.df %>%
    select(matches("^olna_.*")) %>%
    gather(key = "method", value = "n_remove") %>%
    mutate(
      method = str_replace(method, "olna_", ""),
      sim = "Non-parametric",
      lfc_sd = 0.8
    ) ->
  sgdtemp

  bind_rows(pstemp, sgdtemp) %>%
    mutate(sim_lfc = sim) %>%
    mutate(method = case_when(
      grepl("^dout", method) ~ "DESeq2",
      grepl("^eout", method) ~ "edgeR",
      grepl("^nout", method) ~ "NOISeq",
      grepl("^rout", method) ~ "RoPE",
      grepl("^vout", method) ~ "limmav",
      grepl("^wout", method) ~ "Wilcoxon",
    )) %>%
    mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "limmav",
      "NOISeq",
      "Wilcoxon"
    )) -> nremove.dat
  nremove.dat %>%
    group_by(method, sim) %>%
    summarise(m_remove = median(n_remove)) %>%
    ungroup() %>%
    mutate(size = paste0("n/2=", n / 2))
}
