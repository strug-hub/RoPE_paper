here::i_am("code/CF_code/HNE_analysis.R")

library(here)
suppressPackageStartupMessages(library(tidyverse))
library(edgeR)
library(DESeq2)


# Load pre-processed clinical and RNA-Seq counts
load(here("output/CF_data/dat_list_v1.RData"))

Count.v0 <-
  read.table(
    here("output", "CF_data", "active_infection_dge_count.csv"),
    header = T,
    sep = ",",
    row.names = 1
  )

obj.ana3 <- dat.list.v1$gf


y <-
  DGEList(
    counts = counts(obj.ana3),
    group = obj.ana3$longi_pa_binary2
  )


m.formula <-  ~ pc_1 + pc_2 + pc_3 + unc + female + rin + ptprc_tmm + imm_cell_prop +
  longi_pa_binary2

design <- model.matrix(m.formula, data = colData(obj.ana3))

keep <- filterByExpr(y)

y.f <- y[keep,]

count.in <- y.f@.Data[[1]]

# # Filter by expression
# keep <- filterByExpr(y)
# 
# sum(!keep)


# HNE analysis: DESeq2 ----------------------------------------------------
# Design matrix
dds <- DESeqDataSetFromMatrix(countData = count.in,
                              colData = colData(obj.ana3),
                              design = m.formula)

de2res <- DESeq(dds)
de2res.tab <- results(de2res)
de2res.tab.ranked <- de2res.tab[order(de2res.tab$pvalue), ]
# DESeq2 Results
de2res.tab.ranked$Genes <-
  as.character(Count.v0[match(rownames(de2res.tab.ranked), Count.v0$Name), ]$Description)
de2res.tab.ranked.f <- de2res.tab.ranked[!is.na(de2res.tab.ranked$padj),]
res.de2_m3 <- de2res.tab.ranked.f 


res.de2_m3[1:10,]

match("SLC9A3", res.de2_m3$Genes)

sum(res.de2_m3$padj < 0.05)
sum(res.de2_m3$padj < 0.1)
# HNE analysis: edgeR ----------------------------------------------------
y.f <- calcNormFactors(y.f)
y.f <- estimateDisp(y.f, design)
dFit <- glmQLFit(y = y.f , design)
dqlf_0 <- glmQLFTest(dFit)

dqlf <- dqlf_0$table
dqlf$q <- p.adjust(dqlf$PValue, method = "BH")

dqlf <- dqlf[order(dqlf$PValue), ]
dqlf$Genes <-
  as.character(Count.v0[match(rownames(dqlf), Count.v0$Name), ]$Description)
edgeR_m3_qlf <- dqlf

match("SLC9A3", edgeR_m3_qlf$Genes)

edgeR_m3_qlf[1:20,]

edgeR_m3 <- edgeR_m3_qlf
edgeR_m3_qlf[1:10,]
# HNE analysis: voom ----------------------------------------------------
v <- voom(y.f, design, plot = FALSE)

# Test for differential expression
fit <- lmFit(v, design, method = "ls") # default
fit <- eBayes(fit, trend = FALSE, robust = FALSE) # default

# Return top table
de_table <-
  topTable(fit,
           coef = ncol(design),
           n = nrow(y.f),
           sort.by = "P"
  )
de_table$Genes <-
  as.character(Count.v0[match(rownames(de_table), Count.v0$Name), ]$Description)
head(de_table)
lv_m3 <- de_table

match("SLC9A3", lv_m3$Genes)

# HNE analysis: RoPE ----------------------------------------------------
library(roper)

design

RoPE_raw_res <-
  rope(
    datmat = count.in,
    X_model = design,
    x_PI_idx = dim(design)[2]
  )

table_rplr_m3 <- RoPE_raw_res[order(RoPE_raw_res$pvals), ]
table_rplr_m3$Genes <-
  as.character(Count.v0[match(rownames(table_rplr_m3), Count.v0$Name), ]$Description)
RoPE_m3 <- table_rplr_m3

RoPE_m3[1:20,]

match("SLC9A3", RoPE_m3$Genes)

sum(RoPE_m3$padj < 0.05)
sum(RoPE_m3$padj < 0.1)
# save(RoPE_m3, file = here("data_local","out_data","RoPE_HNE.RData"))

RoPE.glist.p05 <- RoPE_m3$Genes[1:sum(RoPE_m3$padj < 0.05)]
RoPE.glist.p1 <- RoPE_m3$Genes[1:sum(RoPE_m3$padj < 0.1)]


# Results in manuscript ---------------------------------------------------
# Tables ------------------------------------------------------------------
slc9a3_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("SLC9A3", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("SLC9A3", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("SLC9A3", lv_m3$Genes), "logFC"],
            RoPE_m3[match("SLC9A3", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("SLC9A3", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("SLC9A3", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("SLC9A3", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("SLC9A3", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("SLC9A3", edgeR_m3$Genes), "q"],
            res.de2_m3[match("SLC9A3", res.de2_m3$Genes), "padj"],
            lv_m3[match("SLC9A3", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("SLC9A3", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("SLC9A3", edgeR_m3$Genes),
    match("SLC9A3", res.de2_m3$Genes),
    match("SLC9A3", lv_m3$Genes),
    match("SLC9A3", RoPE_m3$Genes)
  )
)

knitr::kable(slc9a3_tab, "latex",  digits = c(0, 2, 6,2,2))
slc9a3_tab

formatC(slc9a3_tab$p.val, format = "e", digits = 2)


slc26a4_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("SLC26A4", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("SLC26A4", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("SLC26A4", lv_m3$Genes), "logFC"],
            RoPE_m3[match("SLC26A4", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("SLC26A4", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("SLC26A4", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("SLC26A4", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("SLC26A4", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("SLC26A4", edgeR_m3$Genes), "q"],
            res.de2_m3[match("SLC26A4", res.de2_m3$Genes), "padj"],
            lv_m3[match("SLC26A4", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("SLC26A4", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("SLC26A4", edgeR_m3$Genes),
    match("SLC26A4", res.de2_m3$Genes),
    match("SLC26A4", lv_m3$Genes),
    match("SLC26A4", RoPE_m3$Genes)
  )
)

knitr::kable(slc26a4_tab, "latex",  digits = c(0, 2, 6,2,2))
slc26a4_tab

formatC(slc26a4_tab$p.val, format = "e", digits = 2)


# ANO1

# edgeR_m3[match("ANO1", edgeR_m3$Genes), "logFC"]
ANO1_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("ANO1", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("ANO1", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("ANO1", lv_m3$Genes), "logFC"],
            RoPE_m3[match("ANO1", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("ANO1", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("ANO1", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("ANO1", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("ANO1", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("ANO1", edgeR_m3$Genes), "q"],
            res.de2_m3[match("ANO1", res.de2_m3$Genes), "padj"],
            lv_m3[match("ANO1", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("ANO1", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("ANO1", edgeR_m3$Genes),
    match("ANO1", res.de2_m3$Genes),
    match("ANO1", lv_m3$Genes),
    match("ANO1", RoPE_m3$Genes)
  )
)

knitr::kable(ANO1_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(ANO1_tab, digits=c(0,4,4,4,0))

# APIP
APIP_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("APIP", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("APIP", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("APIP", lv_m3$Genes), "logFC"],
            RoPE_m3[match("APIP", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("APIP", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("APIP", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("APIP", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("APIP", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("APIP", edgeR_m3$Genes), "q"],
            res.de2_m3[match("APIP", res.de2_m3$Genes), "padj"],
            lv_m3[match("APIP", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("APIP", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("APIP", edgeR_m3$Genes),
    match("APIP", res.de2_m3$Genes),
    match("APIP", lv_m3$Genes),
    match("APIP", RoPE_m3$Genes)
  )
)

knitr::kable(APIP_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(APIP_tab, digits=c(0,4,4,4,0))
APIP_tab

# ATP12A
ATP12A_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("ATP12A", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("ATP12A", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("ATP12A", lv_m3$Genes), "logFC"],
            RoPE_m3[match("ATP12A", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("ATP12A", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("ATP12A", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("ATP12A", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("ATP12A", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("ATP12A", edgeR_m3$Genes), "q"],
            res.de2_m3[match("ATP12A", res.de2_m3$Genes), "padj"],
            lv_m3[match("ATP12A", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("ATP12A", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("ATP12A", edgeR_m3$Genes),
    match("ATP12A", res.de2_m3$Genes),
    match("ATP12A", lv_m3$Genes),
    match("ATP12A", RoPE_m3$Genes)
  )
)

knitr::kable(ATP12A_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(ATP12A_tab, digits=c(0,4,4,4,0))

ATP12A_tab


# EHF
EHF_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("EHF", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("EHF", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("EHF", lv_m3$Genes), "logFC"],
            RoPE_m3[match("EHF", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("EHF", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("EHF", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("EHF", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("EHF", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("EHF", edgeR_m3$Genes), "q"],
            res.de2_m3[match("EHF", res.de2_m3$Genes), "padj"],
            lv_m3[match("EHF", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("EHF", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("EHF", edgeR_m3$Genes),
    match("EHF", res.de2_m3$Genes),
    match("EHF", lv_m3$Genes),
    match("EHF", RoPE_m3$Genes)
  )
)

knitr::kable(EHF_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(EHF_tab, digits=c(0,4,4,4,0))
EHF_tab

# MUC4
MUC4_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("MUC4", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("MUC4", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("MUC4", lv_m3$Genes), "logFC"],
            RoPE_m3[match("MUC4", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("MUC4", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("MUC4", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("MUC4", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("MUC4", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("MUC4", edgeR_m3$Genes), "q"],
            res.de2_m3[match("MUC4", res.de2_m3$Genes), "padj"],
            lv_m3[match("MUC4", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("MUC4", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("MUC4", edgeR_m3$Genes),
    match("MUC4", res.de2_m3$Genes),
    match("MUC4", lv_m3$Genes),
    match("MUC4", RoPE_m3$Genes)
  )
)

knitr::kable(MUC4_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(MUC4_tab, digits=c(0,4,4,4,0))
MUC4_tab


# MUC20
MUC20_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("MUC20", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("MUC20", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("MUC20", lv_m3$Genes), "logFC"],
            RoPE_m3[match("MUC20", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("MUC20", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("MUC20", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("MUC20", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("MUC20", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("MUC20", edgeR_m3$Genes), "q"],
            res.de2_m3[match("MUC20", res.de2_m3$Genes), "padj"],
            lv_m3[match("MUC20", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("MUC20", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("MUC20", edgeR_m3$Genes),
    match("MUC20", res.de2_m3$Genes),
    match("MUC20", lv_m3$Genes),
    match("MUC20", RoPE_m3$Genes)
  )
)

knitr::kable(MUC20_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(MUC20_tab, digits=c(0,4,4,4,0))
MUC20_tab

# MUC5AC
MUC5AC_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("MUC5AC", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("MUC5AC", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("MUC5AC", lv_m3$Genes), "logFC"],
            RoPE_m3[match("MUC5AC", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("MUC5AC", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("MUC5AC", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("MUC5AC", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("MUC5AC", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("MUC5AC", edgeR_m3$Genes), "q"],
            res.de2_m3[match("MUC5AC", res.de2_m3$Genes), "padj"],
            lv_m3[match("MUC5AC", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("MUC5AC", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("MUC5AC", edgeR_m3$Genes),
    match("MUC5AC", res.de2_m3$Genes),
    match("MUC5AC", lv_m3$Genes),
    match("MUC5AC", RoPE_m3$Genes)
  )
)

knitr::kable(MUC5AC_tab, "latex",  digits = c(0, 2, 6,2,2))
knitr::kable(MUC20_tab, digits=c(0,4,4,4,0))
MUC5AC_tab

# MUC5B
MUC5B_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("MUC5B", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("MUC5B", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("MUC5B", lv_m3$Genes), "logFC"],
            RoPE_m3[match("MUC5B", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("MUC5B", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("MUC5B", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("MUC5B", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("MUC5B", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("MUC5B", edgeR_m3$Genes), "q"],
            res.de2_m3[match("MUC5B", res.de2_m3$Genes), "padj"],
            lv_m3[match("MUC5B", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("MUC5B", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("MUC5B", edgeR_m3$Genes),
    match("MUC5B", res.de2_m3$Genes),
    match("MUC5B", lv_m3$Genes),
    match("MUC5B", RoPE_m3$Genes)
  )
)

knitr::kable(MUC5B_tab, "latex",  digits = 4)
knitr::kable(MUC5B_tab, digits=c(0,4,4,4,0))

MUC5B_tab


# SLC26A9
SLC26A9_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("SLC26A9", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("SLC26A9", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("SLC26A9", lv_m3$Genes), "logFC"],
            RoPE_m3[match("SLC26A9", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("SLC26A9", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("SLC26A9", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("SLC26A9", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("SLC26A9", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("SLC26A9", edgeR_m3$Genes), "q"],
            res.de2_m3[match("SLC26A9", res.de2_m3$Genes), "padj"],
            lv_m3[match("SLC26A9", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("SLC26A9", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("SLC26A9", edgeR_m3$Genes),
    match("SLC26A9", res.de2_m3$Genes),
    match("SLC26A9", lv_m3$Genes),
    match("SLC26A9", RoPE_m3$Genes)
  )
)

knitr::kable(SLC26A9_tab, "latex",  digits = 4)
knitr::kable(SLC26A9_tab, digits=c(0,4,4,4,0))
SLC26A9_tab


# SLC6A14
SLC6A14_tab <- tibble(
  Method = c("edgeR", "DESeq2", "voom", "RoPE"),
  logFC = c(edgeR_m3[match("SLC6A14", edgeR_m3$Genes), "logFC"],
            log(2 ^ res.de2_m3[match("SLC6A14", res.de2_m3$Genes), "log2FoldChange"]),
            lv_m3[match("SLC6A14", lv_m3$Genes), "logFC"],
            RoPE_m3[match("SLC6A14", RoPE_m3$Genes), "logFC"]),
  p.val = c(edgeR_m3[match("SLC6A14", edgeR_m3$Genes), "PValue"],
            res.de2_m3[match("SLC6A14", res.de2_m3$Genes), "pvalue"],
            lv_m3[match("SLC6A14", lv_m3$Genes), "P.Value"],
            RoPE_m3[match("SLC6A14", RoPE_m3$Genes), "pvals"]),
  adj.p = c(edgeR_m3[match("SLC6A14", edgeR_m3$Genes), "q"],
            res.de2_m3[match("SLC6A14", res.de2_m3$Genes), "padj"],
            lv_m3[match("SLC6A14", lv_m3$Genes), "adj.P.Val"],
            RoPE_m3[match("SLC6A14", RoPE_m3$Genes), "padj"]),
  Rank = c(
    match("SLC6A14", edgeR_m3$Genes),
    match("SLC6A14", res.de2_m3$Genes),
    match("SLC6A14", lv_m3$Genes),
    match("SLC6A14", RoPE_m3$Genes)
  )
)

knitr::kable(SLC6A14_tab, "latex",  digits = 4)
SLC6A14_tab


# RoPE top table ----------------------------------------------------------

RoPE_m3_wr <- RoPE_m3 %>% mutate(edgeR_rank = match(rownames(RoPE_m3), rownames(edgeR_m3)),
                                 DESeq2_rank = match(rownames(RoPE_m3), rownames(res.de2_m3)),
                                 voom_rank = match(rownames(RoPE_m3), rownames(lv_m3)))
RoPE_m3_wr

top_rplr <- RoPE_m3_wr[, c("Genes","logFC","adj","pvals","padj","edgeR_rank","DESeq2_rank","voom_rank")]

top_rplr.p5 <- top_rplr %>% filter(padj < 0.05)

colnames(top_rplr.p5) <- c("Symbol", "logFC", "Adjustment Factor", "p-value", "Adj p-value", "edgeR Rank", "DESeq2 Rank", "voom Rank")

top_rplr.p5$`p-value` <- formatC(top_rplr.p5$`p-value`, format = "e", digits = 2)
formatC(top_rplr.p5$`p-value`, format = "e", digits = 2)

formatC(top_rplr.p5$`Adj p-value`, format = "e", digits = 2)

tt1 <- formatC(top_rplr.p5$`p-value`, format = "e", digits = 2)

tt11 <- gsub("e", " \\\times 10^{", tt1)
tt12 <- gsub("$", "}", tt11)
tt12


# edgeR top table ----------------------------------------------------------

edgeR_m3_top <- edgeR_m3 %>% mutate(RoPE_rank = match(rownames(edgeR_m3), rownames(RoPE_m3)),
                                 DESeq2_rank = match(rownames(edgeR_m3), rownames(res.de2_m3)),
                                 voom_rank = match(rownames(edgeR_m3), rownames(lv_m3)))

edgeR_m3_top[1:20,]


top_edgeR <- edgeR_m3_top[1:20,c("Genes", "logFC", "F","PValue", "q", "RoPE_rank","DESeq2_rank","voom_rank")]

as_tibble(top_edgeR) -> top_edgeR
#top_edgeR$logFC <- num(top_edgeR$logFC, digits = 4)
#top_edgeR$q <- num(top_edgeR$q, digits = 4) 
top_edgeR
formatC(top_edgeR$PValue, format = "f", digits = 4)
top_edgeR$PValue <- formatC(top_edgeR$PValue, format = "f", digits = 4)
top_edgeR$PValue
tt1 <- top_edgeR$PValue
tt11 <- gsub("e", " \\\times 10^{", tt1)
tt12 <- gsub("$", "}", tt11)
tt12
# DESeq2 top table ----------------------------------------------------------
as_tibble(res.de2_m3) -> res.de2_m3.ti

res.de2_m3_top <- res.de2_m3.ti %>% mutate(RoPE_rank = match(rownames(res.de2_m3), rownames(RoPE_m3)),
                                           edgeR_rank = match(rownames(res.de2_m3), rownames(edgeR_m3)),
                                    voom_rank = match(rownames(res.de2_m3), rownames(lv_m3)))

res.de2_m3_top %>% mutate(logFC = log(2^log2FoldChange)) -> res.de2_m3_top

top_DE2 <- res.de2_m3_top[1:20,c("Genes", "logFC", "stat","pvalue", "padj", "RoPE_rank","edgeR_rank","voom_rank")]
#top_DE2$logFC <-  num(top_DE2$logFC, digits = 4)
#top_DE2$padj <- num(top_DE2$padj , digits = 4)
formatC(top_DE2$pvalue, format = "f", digits = 4)
top_DE2$pvalue <- formatC(top_DE2$pvalue, format = "e", digits = 2)
tt1 <- top_DE2$pvalue
tt11 <- gsub("e", " \\\times 10^{", tt1)
tt12 <- gsub("$", "}", tt11)
tt12



# Voom top table ----------------------------------------------------------
lv_m3_top <- lv_m3 %>% mutate(RoPE_rank = match(rownames(lv_m3), rownames(RoPE_m3)),
                              edgeR_rank = match(rownames(lv_m3), rownames(edgeR_m3)),
                                    DESeq2_rank = match(rownames(lv_m3), rownames(res.de2_m3))) %>% as_tibble()
top_lv <- lv_m3_top[1:20,c("Genes", "logFC", "t","P.Value", "adj.P.Val", "RoPE_rank","edgeR_rank","DESeq2_rank")]


#top_lv$logFC <-  num(top_lv$logFC, digits = 4)
#top_lv$adj.P.Val <- num(top_lv$adj.P.Val , digits = 4)
formatC(top_lv$P.Value, format = "f", digits = 4)

top_lv$P.Value <- formatC(top_lv$P.Value, format = "e", digits = 2)

tt1 <- top_lv$P.Value
tt11 <- gsub("e", " \\\times 10^{", tt1)
tt12 <- gsub("$", "}", tt11)
tt12

library(flextable)
ft.er <- flextable(top_edgeR)
ft.er <- colformat_double(x = ft.er,
                      big.mark=",", digits = 4, na_str = "N/A")

ft.de2 <- flextable(top_DE2)
ft.de2 <- colformat_double(x = ft.de2,
                          big.mark=",", digits = 4, na_str = "N/A")

ft.lv <- flextable(top_lv)
ft.lv <- colformat_double(x = ft.lv,
                          big.mark=",", digits = 4, na_str = "N/A")


autofit(ft.er)
autofit(ft.de2)
autofit(ft.lv)

ft.er <- theme_vanilla(ft.er)
ft.de2 <- theme_vanilla(ft.de2)
ft.lv <- theme_vanilla(ft.lv)

save_as_docx("my table 1" = ft.er, "my table 2" = ft.de2,"my table 3" = ft.lv, path = here("output/CF_data/toptabs.docx"))


# Venn diagram of the top 20 ranked DE genes----------------------------------------------------------------
deg_de3 <- res.de2_m3$Genes[1:20]
deg_er3 <- edgeR_m3$Genes[1:20]
deg_lv3 <- lv_m3$Genes[1:20]
deg_rplr3 <- RoPE_m3$Genes[1:20]

Reduce(intersect,list(deg_de3,deg_er3,deg_rplr3))

Reduce(intersect,list(deg_de3,deg_er3,deg_lv3,deg_rplr3))

# library(VennDiagram)
#
# group.venn(
#   list(
#     DESeq2 = deg_de3,
#     edgeR = deg_er3,
#     voom = deg_lv3,
#     RoPE = deg_rplr3
#   ),
#   label = TRUE,
#   cat.pos = c(0, 0, 0, 0),
#   lab.cex = 1.1
# )

library(ggVennDiagram)

gene_list <- list(de = deg_de3, er = deg_er3, lv = deg_lv3, rop = deg_rplr3)

p1 <- ggVennDiagram(gene_list, label = "count", category.names = c("DESeq2","edgeR","voom","RoPE"))
p1

ggsave(file = here("output/old_sim_p/top20.venn.pdf"), p1, width = 6, height = 5)




# GO analysis -------------------------------------------------------------

# out_s <- data.frame(table_rplr_m3$Genes[table_rplr_m3$padj < 0.05])
#
# write.table(
#   out_s,
#   file = here("data_local", "out_data","rope_sig_list.txt"),
#   quote = F,
#   row.names = F,
#   col.names = F
# )

out_s <- data.frame(res.de2_m3$Genes[res.de2_m3$padj < 0.05])

write.table(
  out_s,
  file = here("output/CF_data/DE2_sig_list.txt"),
  quote = F,
  row.names = F,
  col.names = F
)

out_er <- data.frame(top_edgeR$Genes)
out_de2 <- data.frame(top_DE2$Genes)
out_lv <- data.frame(top_lv$Genes)
out_rope <- data.frame(top_rplr.p5$Symbol)

write.table(
  out_er,
  file = here("output/CF_data/er_top_list.txt"),
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  out_de2,
  file = here("output/CF_data/DE2_top_list.txt"),
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  out_lv,
  file = here("output/CF_data/lv_top_list.txt"),
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  out_rope,
  file = here("output/CF_data/rope_top_list.txt"),
  quote = F,
  row.names = F,
  col.names = F
)


RoPE_m3 %>% filter(padj < 0.05) %>%  select(Genes, adj_logLR) %>% mutate(teststat = 2*adj_logLR) %>% select(Genes,teststat) -> gsea.rope

RoPE_m3  %>%  select(Genes, adj_logLR) %>% mutate(teststat = 2*adj_logLR) %>% select(Genes,teststat) -> gsea.rope.full


write.table(
  gsea.rope,
  file = here("output/CF_data/rope_top_gesa.txt"),
  sep = '\t',
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  gsea.rope.full,
  file = here("output/CF_data/rope_top_gesa_full.txt"),
  sep = '\t',
  quote = F,
  row.names = F,
  col.names = F
)


# GSEA results ------------------------------------------------------------
GSEA_res <- read.table(here("output/CF_data/gsea_report_for_na_pos_1678145357326.tsv"),sep = '\t',header = T) %>%  as_tibble()
colnames(GSEA_res) 

GSEA_res %>% slice_head(n=10) %>% select(NAME, SIZE, ES,NES,NOM.p.val ,FDR.q.val) -> GSEA_res.out

GSEA_res.ft <- flextable(GSEA_res.out)

GSEA_res.ft <- colformat_double(x = GSEA_res.ft,
                          big.mark=",", digits = 4, na_str = "N/A")

GSEA_res.ft <- theme_vanilla(GSEA_res.ft)

save_as_docx("my table 1" = ft.er, "my table 2" = ft.de2,"my table 3" = ft.lv, "my table 4" = GSEA_res.ft,path = here("output/CF_data/toptabs.docx"))
