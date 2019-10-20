library(DESeq2)
library("BiocParallel")
register(MulticoreParam(4))
Matrix_Total = read.csv('tumoral_healthy.txt', sep = '\t', header = T)
rownames(Matrix_Total) <- Matrix_Total$Gene
Matrix_Total$Gene <- NULL

Matrix_Total$TCGA.A6.6650.01A_x <- NULL #1
Matrix_Total$TCGA.A6.6650.01A_y <- NULL
Matrix_Total$TCGA.A6.6650.01B <- NULL 
Matrix_Total$TCGA.A6.5659.01A_x <- NULL #2
Matrix_Total$TCGA.A6.5659.01A_y <- NULL
Matrix_Total$TCGA.A6.5659.01B <- NULL
Matrix_Total$TCGA.A6.2684.01A_x <- NULL #3
Matrix_Total$TCGA.A6.2684.01A_y <- NULL 
Matrix_Total$TCGA.A6.2684.01C <- NULL
Matrix_Total$TCGA.A6.5656.01A_x <- NULL #4
Matrix_Total$TCGA.A6.5656.01A_y <- NULL
Matrix_Total$TCGA.A6.5656.01B <- NULL
Matrix_Total$TCGA.A6.3810.01A_x <- NULL #5
Matrix_Total$TCGA.A6.3810.01A_y <- NULL
Matrix_Total$TCGA.A6.3810.01B <- NULL
Matrix_Total$TCGA.A6.2672.01A <- NULL #6
Matrix_Total$TCGA.A6.2672.01B <- NULL
Matrix_Total$TCGA.A6.5661.01A <- NULL #7
Matrix_Total$TCGA.A6.5661.01B <- NULL
Matrix_Total$TCGA.A6.5665.01A <- NULL #8
Matrix_Total$TCGA.A6.5665.01B <- NULL
Matrix_Total$TCGA.A6.2677.01A <- NULL #9
Matrix_Total$TCGA.A6.2677.01B <- NULL

iqrs = apply(Matrix_Total, 1, IQR)	# by setting the second argument to 1, we are specifying to run IQR by row in the matrix. How would you modify it to run the IQR function by columns?
th = median(iqrs)
idx = iqrs >= th
filt_exprdata = Matrix_Total[idx,]  




vector_t_test = c()   ##1 tumoral 0 healthy
library(stringr)
values_tumoral = str_pad(1:9, pad = 0,width = 2 , "left")
values_healthy = str_pad(10:19, pad = 0,width = 2 , "left")
ids = substr(colnames(filt_exprdata), start = 14, stop = 15)
for (col in ids ){
  if (col%in% values_tumoral){
    vector_t_test = c(vector_t_test, 'Treated')
  }else{
    vector_t_test = c(vector_t_test, 'Untreated')
  }
}
st = as.factor(vector_t_test)
coldata <- data.frame(row.names=colnames(filt_exprdata), st)
ddshTSeq <- DESeqDataSetFromMatrix(countData =filt_exprdata, design = ~ st, colData = coldata)
keep <- rowSums(counts(ddshTSeq)) >= 10
dds <- ddshTSeq[keep,]
dds <- DESeq(ddshTSeq, parallel = TRUE)
dds_results <- results(dds)
dds_results_df <- as.data.frame(dds_results)
write.table(dds_results_df, 'DESeq2_IQR_Median_results.csv', quote = F, sep = '\t', row.names = T, col.names = T)


library(EnhancedVolcano)

EnhancedVolcano(dds_results,
                lab = rownames(dds_results),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5)

dds_results_filtered = dds_results_df[dds_results_df$padj<0.05,] ## 26962 genes
## if we filter for the log2fold change we get 4235 genes is we filter for |log2fc| > 2
## we show here the volcano plot.

EnhancedVolcano(dds_results_filtered_lg2fc,
                lab = rownames(dds_results_filtered_lg2fc),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=colnames(coldata)) 

