### imports

library(edgeR)
library(lattice)
library(biomaRt)
library(gplots)
library(gridExtra)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(Glimma)

load("course.Rdata")

groups <- factor(c("PBS","EGF","PBS","EGF"))

label <- c("PBS2", "EGF1", "PBS1", "EGF2")

print(colnames(subread_counts$counts))
print(colnames(subread_counts$counts))

DGEObj <- DGEList(group=groups, counts=subread_counts$counts, genes=subread_counts$annotation[,c("GeneID","Length")])

print(DGEObj)

design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
print(design)

contrasts <- makeContrasts(
  EGF_vs_PBS = EGF - PBS,
  levels=design
)
counts <- DGEObj$counts

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

fit <- glmFit(DGEObj, design) # fit generalized linear model
lrt <- glmLRT(fit, contrast = contrasts[, "EGF_vs_PBS"])

# main data table of edgeR results
edgerTable <- topTags(lrt, n = nrow(DGEObj))$table

####################################################################################
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

ttop_dge <- edgerTable
head(ttop_dge)