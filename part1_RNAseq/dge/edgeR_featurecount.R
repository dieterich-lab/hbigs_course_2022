#!/usr/bin/env Rscript

### imports

library(edgeR)
library(lattice)
library(biomaRt)
library(gplots)
library(gridExtra)
library(openxlsx)
library(tools)
library(ggplot2)
library(dplyr)
library(Glimma)

# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

contrast <- args[1]
global_pval <- as.numeric(args[2])
global_fdr <- as.numeric(args[3])
log_threshold <- as.numeric(args[4])
threshold <- as.numeric(args[5])
focus <- args[6]


## fix "not find zip" error from openxlsx
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
require(openxlsx)

# visible R code
options(echo = TRUE)

# set base output directory
baseDir <- paste("./edgeR_p", global_pval, "__FDR_", global_fdr, "__thresh_", threshold, "/", sep = "")
dir.create(baseDir, showWarnings = TRUE, recursive = TRUE, mode = "0755")



load("course.Rdata")


############
source(file = paste(contrast, ".R", sep = ""))
###########

head(DGEObj$samples)

rownames(DGEObj$samples) <- gsub("\\.", "_", rownames(DGEObj$samples))
rownames(DGEObj$samples) <- gsub("KD_", "", rownames(DGEObj$samples))
rownames(DGEObj$samples) <- gsub("_circSLC8A1", "", rownames(DGEObj$samples))
rownames(DGEObj$samples) <- gsub("_bam", "", rownames(DGEObj$samples))
rownames(DGEObj$counts) <- gsub("_1_mRNA_Seq_PBS", "", rownames(DGEObj$counts))
rownames(DGEObj$counts) <- gsub("_1_mRNA_Seq_EGF", "", rownames(DGEObj$counts))
rownames(DGEObj$counts) <- gsub("_2_mRNA_Seq_PBS", "", rownames(DGEObj$counts))
rownames(DGEObj$counts) <- gsub("_2_mRNA_Seq_EGF", "", rownames(DGEObj$counts))

colnames(DGEObj$counts) <- gsub("\\.", "_", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("KD_", "", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("_1_mRNA_Seq_PBS", "", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("_1_mRNA_Seq_EGF", "", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("_2_mRNA_Seq_PBS", "", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("_2_mRNA_Seq_EGF", "", colnames(DGEObj$counts))
colnames(DGEObj$counts) <- gsub("_bam", "", colnames(DGEObj$counts))


DGEObj$samples
head(DGEObj$counts)
head(DGEObj$genes)
DGEObj



# keep only candidates with > 2 counts per million mapped reads AND in 3 libs
keep <- rowSums(cpm(DGEObj) > 1) >= 3

# select sub set from above
DGEObj <- DGEObj[keep, keep.lib.sizes = FALSE]

dim(DGEObj)

DGEObj <- calcNormFactors(DGEObj) # recalc norm factors


DGEObj <- estimateDisp(DGEObj, design) # estimate dispersion
DGEObj$common.dispersion

# MDS plot
pdf(paste(baseDir, "plot_mds.pdf", sep = ""), title = "MA plot")
col.status <- c("blue", "red", "green", "cyan", "black")[DGEObj$samples$group]
plotMDS(DGEObj, col = col.status)

dev.off()

fit <- glmFit(DGEObj, design) # fit generalized linear model
lrt <- glmLRT(fit, contrast = contrasts[, contrast])

isDE <- as.logical(decideTestsDGE(lrt, p.value = global_pval))
DEnames <- rownames(DGEObj)[isDE]
summary(isDE)

# main data table of edgeR results
edgerTable <- topTags(lrt, n = nrow(DGEObj))$table




y <- cpm(DGEObj, log = TRUE, prior.count = 1)
lcpm <- cpm(DGEObj, log = TRUE)


glMDSPlot(lcpm, labels = paste(groups, sep = "_"), groups = DGEObj$samples,
  launch = FALSE, path = baseDir)


####################################################################################
pdf(paste(baseDir, contrast, ".pdf", sep = ""), title = paste(contrast, " results"))
plotSmear(lrt, xlab = "Log Concentration" , ylab = "Log Fold-Change", smooth.scatter = FALSE, lowess = FALSE, de.tags = DEnames, pch = 19, cex = 0.4, main = paste(contrast, " expression", sep = ""))
abline(h = c((- 1) * threshold, threshold), col = "blue")

dev.off()



edgerResults <- edgerTable[edgerTable$PValue <= global_pval ,]


## done excel export
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="ensembl.org")

genemap <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), values = edgerTable$GeneID, mart = ensembl)
genemap2 <- getBM(attributes = c("ensembl_gene_id", "go_id"), filter = "ensembl_gene_id", values = edgerTable$GeneID, mart = ensembl)

edgerTable.idx <- match(edgerTable$GeneID, genemap$ensembl_gene_id)
edgerTable.idx2 <- match(edgerTable$GeneID, genemap2$ensembl_gene_id)

# add biomart data to results tablesTRUE
edgerTable$external_gene_name <- genemap$external_gene_name[edgerTable.idx]
# edgerTable$entrezgene <- genemap$entrezgene[edgerTable.idx]
edgerTable$description <- genemap$description[edgerTable.idx]
edgerTable$GO <- genemap2$go_id[edgerTable.idx2]
edgerTable$ensembl_gene_id <- genemap$ensembl_gene_id[edgerTable.idx]
lrt$table$gene <- genemap$external_gene_name[edgerTable.idx]
gene_names <- data.frame(lrt$table$gene)

edgerTable <- edgerTable[c("external_gene_name", "GeneID" , "Length" , "logFC", "logCPM", "PValue", "FDR", "GO", "description")]


edgerResults <- edgerTable[edgerTable$PValue < global_pval,]
edgerFiltered <- edgerResults[edgerResults$FDR <= global_fdr & abs(edgerResults$logFC) > log_threshold ,]
edgerSortedByFC <- edgerResults[order(edgerFiltered$logFC),]
edgerSortedByFCselection <- edgerSortedByFC[edgerSortedByFC$logFC <= (- 1) * log_threshold | edgerSortedByFC$logFC >= log_threshold,]

glMDPlot(lrt, anno=gene_names, status = decideTestsDGE(lrt, p.value = global_pval),
  counts = DGEObj, groups = groups, transform = TRUE, launch = FALSE, path = baseDir,
  side.main = "Symbol", main = "Gene expression")


## create excel workbook and file
wb <- createWorkbook()

addWorksheet(wb, sheetName = paste("p<", global_pval, sep = ""));
writeDataTable(wb, sheet = 1, x = edgerResults, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR<", global_fdr, sep = ""));
writeDataTable(wb, sheet = 2, x = edgerFiltered, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR<", global_fdr, ", by FC", sep = ""));
writeDataTable(wb, sheet = 3, x = edgerSortedByFC, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR< ", global_fdr, ", logFC>", log_threshold, sep = ""));
writeDataTable(wb, sheet = 4, x = edgerSortedByFCselection, rowNames = FALSE);

saveWorkbook(wb, paste(output = baseDir, contrast, "_edger.xlsx", sep = ""), overwrite = TRUE)





