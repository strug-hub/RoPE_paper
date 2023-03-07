here::i_am("code/figure_code/draw_benmark_time.R")

library(here)

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggsci)
library(patchwork)

benmark.n50 <- readRDS(here("output/bench_time/nc.1_pn0.9_ns.50_it.3_fil.FALSE.btime.new.RDS"))
benmark.n100 <- readRDS(here("output/bench_time/nc.1_pn0.9_ns.100_it.3_fil.FALSE.btime.new.RDS"))
benmark.n200 <- readRDS(here("output/bench_time/nc.1_pn0.9_ns.200_it.3_fil.FALSE.btime.new.RDS"))
benmark.n400 <- readRDS(here("output/bench_time/nc.1_pn0.9_ns.400_it.3_fil.FALSE.btime.new.RDS"))

get.time <- function(ben.list, n) {
  tt <- sapply(ben.list, FUN = function(obj) {
    unlist(obj$median)
  })
  tt <- as_tibble(tt) %>%
    mutate(method = c("limmav", "DESeq2", "edgeR", "RoPE", "Wilcoxon", "NOISeq", "dearseq")) %>%
    mutate(n = paste0("n/2=", n / 2))
  return(tt)
}

time.dat <- bind_rows(
  get.time(benmark.n50, n = 50),
  get.time(benmark.n100, n = 100),
  get.time(benmark.n200, n = 200),
  get.time(benmark.n400, n = 400)
) %>%
  mutate(method = fct_relevel(
    method,
    "RoPE",
    "edgeR",
    "DESeq2",
    "limmav",
    "NOISeq",
    "Wilcoxon",
    "dearseq"
  )) %>%
  mutate(n = fct_relevel(n, "n/2=25", "n/2=50", "n/2=100", "n/2=200"))

p0 <- ggplot(time.dat, aes(x = method, y = value)) +
  geom_col() +
  facet_grid(. ~ n) +
  theme_bw() +
  ylab("Computational Time (seconds)")

ggsave(file = here("output/figures_time/time.0.pdf"), p0, width = 15, height = 8)

bold.italic.text <- element_text(face = "bold.italic")

p1 <- ggplot(time.dat %>% filter(n == "n/2=25"), aes(x = method, y = value)) +
  ggtitle("n/2=25") + 
  scale_y_continuous(limits = c(0, 550), breaks = seq(0,500,100)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p2 <- ggplot(time.dat %>% filter(n == "n/2=50"), aes(x = method, y = value)) +
  ggtitle("n/2=50") + 
  scale_y_continuous(limits = c(0, 550), breaks = seq(0,500,100)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p3 <- ggplot(time.dat %>% filter(n == "n/2=100"), aes(x = method, y = value)) +
  ggtitle("n/2=100") + 
  scale_y_continuous(limits = c(0, 550), breaks = seq(0,500,100)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p4 <- ggplot(time.dat %>% filter(n == "n/2=200"), aes(x = method, y = value)) +
  ggtitle("n/2=200") + 
  scale_y_continuous(limits = c(0, 550), breaks = seq(0,500,100)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)



p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "Computational Time (seconds)")

ggsave(file = here("output/figures_time/time.all.pdf"), g, width = 18, height = 6)




time.dat.1 <- time.dat %>% filter(method != "dearseq") %>% droplevels()


p1 <- ggplot(time.dat.1 %>% filter(n == "n/2=25"), aes(x = method, y = value)) +
  ggtitle("n/2=25") + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0,80,20)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p2 <- ggplot(time.dat.1 %>% filter(n == "n/2=50"), aes(x = method, y = value)) +
  ggtitle("n/2=50") + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0,80,20)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p3 <- ggplot(time.dat.1 %>% filter(n == "n/2=100"), aes(x = method, y = value)) +
  ggtitle("n/2=100") + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0,80,20)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)
p4 <- ggplot(time.dat.1 %>% filter(n == "n/2=200"), aes(x = method, y = value)) +
  ggtitle("n/2=200") + 
  scale_y_continuous(limits = c(0, 80), breaks = seq(0,80,20)) + 
  geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = bold.italic.text)



p.res <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")

gt <- patchwork::patchworkGrob(p.res)
g <- gridExtra::grid.arrange(gt, left = "Computational Time (seconds)")


ggsave(file = here("output/figures_time/time.ndear.pdf"), g, width = 18, height = 6)

