
doDGE <- function(rawct, annot){
  dds <- DESeqDataSetFromMatrix(countData = rawct, colData = annot, design = ~ Purist)
  dds <- DESeq(dds)
  res <- results(dds)
  significant_genes <- res$stat[which(res$padj < 0.05)]
  ensg = rownames(res)[which(res$padj < 0.05)]
  names(significant_genes) = ensg
  # names(significant_genes) = geneannot$GeneName[match(ensg, geneannot$GeneID)]
  return(res)
}
