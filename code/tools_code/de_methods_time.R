##################
## DE functions
##################
## design_mat should have the intercept term in all of these

get_voom <- function(countdat, design_mat) {
  vout <- limma::voom(counts = countdat, design = design_mat)
  lout <- limma::lmFit(vout)
  eout <- limma::eBayes(lout)
  retlist <- list(
    bhat = coefficients(eout)[, 2],
    pval = eout$p.value[, 2],
    qval = p.adjust(eout$p.value[, 2], method = "BH")
    # n.remove = sum(is.na(eout$p.value))
  )
  return(retlist)
}

get_DESeq2 <- function(countdat, design_mat) {
  colnames(design_mat) <- c("I", "D")
  design_mat <- as.data.frame(design_mat)
  design_mat$D <- factor(design_mat$D)
  trash <- capture.output({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdat, colData = design_mat, design = ~D)
    dds <- DESeq2::DESeq(object = dds, quiet = TRUE)
    res <- DESeq2::results(dds)
  })
  retlist <- list(bhat = res$log2FoldChange, 
                  pval = res$pvalue, 
                  qval = res$padj
                  # n.remove = sum(is.na(res$padj))
                  )
  return(retlist)
}


# get_edgeR <- function(countdat, design_mat) {
#   trash <- capture.output({
#     dge <- edgeR::DGEList(counts = countdat, group = design_mat[, 2])
#     dge <- edgeR::estimateDisp(y = dge)
#     efit <- edgeR::exactTest(dge)
#   })
#   retlist <- list(bhat = efit$table$logFC,
#                   pval = efit$table$PValue)
# }


get_rope <- function(countdat, design_mat) {
  trash <- capture.output({
    res <- roper::rope(datmat = countdat, X_model = design_mat, x_PI_idx = dim(design_mat)[2])
  })
  retlist <- list(
    bhat = log2(exp(res$logFC)),
    pval = res$pvals,
    qval = res$padj
    # n.remove = sum(is.na(res$padj))
  )
  return(retlist)
}



get_edgeR <- function(countdat, design_mat) {
  trash <- capture.output({
    dge <- edgeR::DGEList(counts = countdat, group = design_mat[, 2])
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(y = dge, design = design_mat)
    fit <- edgeR::glmQLFit(dge)
    efit <- edgeR::glmQLFTest(fit)
  })
  retlist <- list(
    bhat = efit$table$logFC,
    pval = efit$table$PValue,
    qval = p.adjust(efit$table$PValue, method = "BH")
    # n.remove = sum(is.na(efit$table$PValue))
  )
  return(retlist)
}

get_wilcoxon <- function(countdat, design_mat) {
  trash <- capture.output({
    dge <- edgeR::DGEList(counts = countdat, group = design_mat[, 2])
    # keep <- filterByExpr(y)
    # y <- y[keep,keep.lib.sizes=FALSE]
    dge <- edgeR::calcNormFactors(dge)
    count_norm <- edgeR::cpm(dge)
    count_norm <- as.data.frame(count_norm)

    dataMem1 <- count_norm[, design_mat[, 2] == 0]
    dataMem2 <- count_norm[, design_mat[, 2] == 1]
    suppressWarnings(pvalue <- matrixTests::row_wilcoxon_twosample(dataMem1, dataMem2)$pvalue)
  })
  retlist <- list(
    bhat = rep(NA, nrow(countdat)),
    pval = pvalue,
    qval = p.adjust(pvalue, method = "BH")
    # n.remove = sum(is.na(pvalue))
  )
  return(retlist)
}

get_NOISeq <- function(countdat, design_mat) {
  trash <- capture.output({
    conditions <- factor(design_mat[, 2])
    data <- NOISeq::readData(data = countdat, factors = as.data.frame(conditions))
    res <- NOISeq::noiseqbio(data, k = 0.5, norm = "tmm", factor = "conditions", random.seed = 12345, filter = 0, cv.cutoff = 100, cpm = 1)
  })
  retlist <- list(
    bhat = res@results[[1]]$log2FC,
    pval = rep(NA, nrow(countdat)),
    qval = 1 - res@results[[1]]$prob
    # n.remove = sum(is.na(res@results[[1]]$prob))
  )
  return(retlist)
}

get_dearseq <- function(countdat, design_mat){
  trash <- capture.output({
    conditions<-matrix(as.numeric(design_mat[,2]),ncol=1)
    res=
      suppressMessages(suppressWarnings(dearseq::dear_seq(exprmat=as.matrix(countdat), variables2test=conditions, 
                          parallel_comp=F, preprocessed=F, na.rm_dearseq = T)))
  })
  retlist <- list(
    bhat = rep(NA, nrow(res$pvals)),
    pval = res$pvals$rawPval,
    qval = res$pvals$adjPval
    # n.remove = nrow(countdat) - nrow(res$pvals)
  )
}
