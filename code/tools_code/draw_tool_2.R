here::i_am("code/tools_code/draw_tool_2.R")

library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggsci)


# fdp_plots ---------------------------------------------------------------
draw_fdp_a <- function(NP.df, NB.df, q = 0.1, n = 50) {
  NB.df %>%
    select(matches("^fpr_.*"), lfc_sd) %>%
    select(ends_with(paste0("q_", q)), lfc_sd) %>%
    gather(-lfc_sd, key = "method", value = "fdp") %>%
    mutate(
      method = str_replace(method, "fpr_", ""),
      sim = "NB-parametric"
    ) ->
    pstemp
  
  NP.df %>%
    select(matches("^fpr_.*")) %>%
    select(ends_with(paste0("q_", q))) %>%
    gather(key = "method", value = "fdp") %>%
    mutate(
      method = str_replace(method, "fpr_", ""),
      sim = "Non-parametric",
      lfc_sd = 0.8
    ) ->
    sgdtemp
  
  bind_rows(pstemp, sgdtemp) %>%
    mutate(sim_lfc = sim) %>%
    mutate(method = case_when(
      grepl("^dout", method) ~ "DESeq2",
      grepl("^eout", method) ~ "edgeR",
      grepl("^nout", method) ~ "NOISeq_nf",
      grepl("^rout", method) ~ "RoPE",
      grepl("^vout", method) ~ "limmav",
      grepl("^wout", method) ~ "Wilcoxon",
      grepl("^dnfout", method) ~ "DESeq2_nf",
      grepl("^nfout", method) ~ "NOISeq",
      grepl("^xout", method) ~ "dearseq"
    )) %>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    ))  %>%
    ggplot(aes(x = method, y = fdp, color = sim_lfc)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_color_colorblind(name = "Simulation\nMethod") +
    theme_bw() +
    geom_hline(yintercept = q, lty = 2) +
    xlab("Method") +
    ylab("FDP") +
    ggtitle(paste0("n/2=", n / 2, ", q=", q)) ->
    pl
  pl
  return(pl)
}
# power plots -------------------------------------------------------------
draw_power_a <- function(NP.df, NB.df, q = 0.1, n = 50) {
  NB.df %>%
    select(matches("^power_.*"), lfc_sd) %>%
    select(ends_with(paste0("q_", q)), lfc_sd) %>%
    gather(-lfc_sd, key = "method", value = "power") %>%
    mutate(
      method = str_replace(method, "power_", ""),
      sim = "NB-parametric"
    ) ->
    pstemp
  
  NP.df %>%
    select(matches("^power_.*")) %>%
    select(ends_with(paste0("q_", q))) %>%
    gather(key = "method", value = "power") %>%
    mutate(
      method = str_replace(method, "power_", ""),
      sim = "Non-parametric",
      lfc_sd = 0.8
    ) ->
    sgdtemp
  
  bind_rows(pstemp, sgdtemp) %>%
    mutate(sim_lfc = sim) %>%
    mutate(method = case_when(
      grepl("^dout", method) ~ "DESeq2",
      grepl("^eout", method) ~ "edgeR",
      grepl("^nout", method) ~ "NOISeq_nf",
      grepl("^rout", method) ~ "RoPE",
      grepl("^vout", method) ~ "limmav",
      grepl("^wout", method) ~ "Wilcoxon",
      grepl("^dnfout", method) ~ "DESeq2_nf",
      grepl("^nfout", method) ~ "NOISeq",
      grepl("^xout", method) ~ "dearseq"
    )) %>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    ))%>%
    ggplot(aes(x = method, y = power, color = sim_lfc)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_color_colorblind(name = "Simulation\nMethod") +
    theme_bw() +
    xlab("Method") +
    ylab("power") +
    ggtitle(paste0("n/2=", n / 2, ", q=", q)) ->
    pl
  return(pl)
}


# grid --------------------------------------------------------------------

# draw_power_grid --------------------------------------------------------------
draw_power_grid_a <- function(NP.df, NB.df, n){
  
  NB.df %>%
    select(matches("^power_.*")) %>%
    mutate(iter = 1:nrow(NB.df)) %>%
    select(matches("^power_.*"), iter) %>%
    pivot_longer(cols = !iter, names_to = "column_name", values_to = "power") %>%
    separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
    mutate(method = gsub("power_", "", method)) %>%
    mutate(method = factor(recode(method,
                                  "dout" = "DESeq2",
                                  "eout" = "edgeR",
                                  "vout" = "limmav",
                                  "rout" = "RoPE",
                                  "wout" = "Wilcoxon",
                                  "nout" = "NOISeq_nf",
                                  "dnfout" = "DESeq2_nf",
                                  "nfout" = "NOISeq",
                                  "xout" = "dearseq"
                                  
    )))%>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    )) %>%
    group_by(method, qvalue) %>%
    summarise(med_power = median(power)) %>%
    ungroup() %>%
    mutate(qvalue = as.numeric(qvalue)) -> power.nb
  
  
  NP.df %>%
    select(matches("^power_.*")) %>%
    mutate(iter = 1:nrow(NB.df)) %>%
    select(matches("^power_.*"), iter) %>%
    pivot_longer(cols = !iter, names_to = "column_name", values_to = "power") %>%
    separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
    mutate(method = gsub("power_", "", method)) %>%
    mutate(method = factor(recode(method,
                                  "dout" = "DESeq2",
                                  "eout" = "edgeR",
                                  "vout" = "limmav",
                                  "rout" = "RoPE",
                                  "wout" = "Wilcoxon",
                                  "nout" = "NOISeq_nf",
                                  "dnfout" = "DESeq2_nf",
                                  "nfout" = "NOISeq",
                                  "xout" = "dearseq"
                                  
    )))%>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    )) %>%
    group_by(method, qvalue) %>%
    summarise(med_power = median(power)) %>%
    ungroup() %>%
    mutate(qvalue = as.numeric(qvalue)) -> power.np
  
  power.dat <- bind_rows(
    power.nb %>% mutate(sim = "NB-parametric"),
    power.np %>% mutate(sim = "Non-parametric")
  )
  
  pl <- ggplot(power.dat, aes(x = qvalue, y = med_power, color = method)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    facet_grid(sim ~ .) + ggtitle(paste0("n/2=", n / 2))
  
  pl
  return(pl)
}

# draw_fpr_grid -----------------------------------------------------------
draw_fpr_grid_a <- function(NP.df, NB.df, n){
  NB.df %>%
    select(matches("^fpr_.*")) %>% mutate(iter = 1:nrow(NB.df)) %>%
    select(matches("^fpr_.*"),iter) %>%
    pivot_longer(cols = !iter, names_to = "column_name", values_to = "fpr") %>%
    separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
    mutate(method = gsub("fpr_", "", method)) %>%
    mutate(method = factor(recode(method,
                                  "dout" = "DESeq2",
                                  "eout" = "edgeR",
                                  "vout" = "limmav",
                                  "rout" = "RoPE",
                                  "wout" = "Wilcoxon",
                                  "nout" = "NOISeq_nf",
                                  "dnfout" = "DESeq2_nf",
                                  "nfout" = "NOISeq",
                                  "xout" = "dearseq"
                                  
    )))%>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    )) %>% group_by(method, qvalue) %>% summarise(med_fpr = median(fpr)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue)) -> fpr.nb
  
  
  
  NP.df %>%
    select(matches("^fpr_.*")) %>% mutate(iter = 1:nrow(NB.df)) %>%
    select(matches("^fpr_.*"),iter) %>%
    pivot_longer(cols = !iter, names_to = "column_name", values_to = "fpr") %>%
    separate(column_name, into = c("method", "qvalue"), sep = "_q_") %>%
    mutate(method = gsub("fpr_", "", method)) %>%
    mutate(method = factor(recode(method,
                                  "dout" = "DESeq2",
                                  "eout" = "edgeR",
                                  "vout" = "limmav",
                                  "rout" = "RoPE",
                                  "wout" = "Wilcoxon",
                                  "nout" = "NOISeq_nf",
                                  "dnfout" = "DESeq2_nf",
                                  "nfout" = "NOISeq",
                                  "xout" = "dearseq"
                                  
    )))%>% mutate(method = fct_relevel(
      method,
      "RoPE",
      "edgeR",
      "DESeq2",
      "DESeq2_nf",
      "limmav",
      "NOISeq",
      "NOISeq_nf",
      "Wilcoxon",
      "dearseq"
    )) %>% group_by(method, qvalue) %>% summarise(med_fpr = median(fpr)) %>% ungroup() %>% mutate(qvalue = as.numeric(qvalue)) -> fpr.np
  
  
  
  fpr.dat <- bind_rows(fpr.nb %>% mutate(sim = "NB-parametric"),
                       fpr.np %>% mutate(sim = "Non-parametric"))
  
  pl <- ggplot(fpr.dat, aes(x = qvalue, y = med_fpr, color = method)) + geom_line(linewidth = 0.8) + theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    facet_grid(sim~.)  + ggtitle(paste0("n/2=", n / 2))
  
  return(pl)
}




