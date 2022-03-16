set.seed(1234)

library("readr")
library("vsn")
library("dplyr")
library("limma")
library("ggplot2")
library("ggrepel")
library("BioNet")
library("igraph")
library("OmnipathR")
library("ggpubr")
library("mixOmics")
library("M2SMF")
library("SNFtool")
library("NEMO")
library("fgsea")
library("GSA")
library("VennDiagram")
library("RColorBrewer")
library("ggVennDiagram")

# Differential Gene Expression Data
load(file = "../Data/ttop_dge.RData")
head(ttop_dge[, 1:(ncol(ttop_dge)-1)])

# Differential Protein Abundance
load(file = "../Data/ttop_prot.RData")
head(ttop_prot)

# Processed Gene Expression Data across EGF and PBS samples at time-point 60min
load(file = "../Data/proc_gene_data.RData")
head(proc_gene_data)

# Processed Protein Abundance Data across EGF and PBS samples at time-point 60min
load(file = "../Data/proc_prot_data.RData")
head(proc_prot_data)

# # Order by FDR and filter duplicates
# ttop_dge <- ttop_dge[order(ttop_dge$FDR, decreasing = FALSE), ]
# ttop_prot <- ttop_prot[order(ttop_prot$EGF_60_vs_PBS_60_p.adj, decreasing = FALSE), ]
# 
# ttop_dge <- ttop_dge[-which(duplicated(ttop_dge$external_gene_name)), ]
# ttop_prot <- ttop_prot[-which(duplicated(ttop_prot$Gene)), ]

# We find common Genes and and filter each data
common_genes <- intersect(x = ttop_dge$external_gene_name, y = ttop_prot$Gene)
dge <- ttop_dge[which(ttop_dge$external_gene_name%in%common_genes), ]
prot <- ttop_prot[which(ttop_prot$Gene%in%common_genes), ]

# We create the data-frame for plotting the correlation
data <- matrix(data = , nrow = length(common_genes), ncol = 2)
rownames(data) <- common_genes[order(common_genes)]
colnames(data) <- c("diff_genes", "diff_prot")
data[, 1] <- dge$logFC[order(dge$external_gene_name)]
data[, 2] <- prot$EGF_60_vs_PBS_60_diff[order(prot$Gene)]

data <- as.data.frame(data)
head(data)

# We do ascatter plot of gene expression and protein abundance and estimate the 
# Pearson correlation between them
sp <- ggscatter(data, x = "diff_genes", y = "diff_prot", #mention data and axis 
                add = "reg.line",  # Add regression line
                add.params = list(color = "red", fill = "lightgray"), # Customize regression line
                conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 3, label.y = 30)# Add correlation coefficient
sp

## Correlation Plot
cim(cor(t(proc_gene_data), 
        t(proc_prot_data)), 
    xlab = "proteins", ylab = "genes")

# Filter for Common Genes
common_genes <- intersect(x = ttop_dge$external_gene_name, y = ttop_prot$Gene)
omic1 <- ttop_dge[which(ttop_dge$external_gene_name%in%common_genes), ]
omic1 <- as.matrix(omic1$logFC[order(omic1$external_gene_name)])
omic2 <- ttop_prot[which(ttop_prot$Gene%in%common_genes), ]
omic2 <- as.matrix(omic2$EGF_60_vs_PBS_60_diff[order(omic2$Gene)])
rownames(omic1) <- common_genes[order(common_genes)]
rownames(omic2) <- common_genes[order(common_genes)]

omics.list = list(t(omic1), t(omic2))

# Supervised Clustering - 2 Clusters
clustering = nemo.clustering(omics.list = omics.list, num.clusters = 2)

data <- cbind(omic1, omic2, as.matrix(clustering))
colnames(data) <- c("diff_genes", "diff_prot", "cluster")
data <- as.data.frame(data)
data$cluster <- as.factor(data$cluster)
head(data)

p <- ggplot(data, aes(x=diff_genes, y=diff_prot, color=cluster, shape=cluster)) +
  geom_point() + 
  theme_classic() +
  ggtitle("Scatter Plot and Clustering")
p+scale_color_brewer(palette="Dark2")

# Unsupervised Clustering
clustering = nemo.clustering(omics.list = omics.list, num.clusters = NA)

data <- cbind(omic1, omic2, as.matrix(clustering))
colnames(data) <- c("diff_genes", "diff_prot", "cluster")
data <- as.data.frame(data)
data$cluster <- as.factor(data$cluster)
head(data)

p <- ggplot(data, aes(x=diff_genes, y=diff_prot, color=cluster, shape=cluster)) +
  geom_point() + 
  theme_classic() +
  ggtitle("Scatter Plot and Clustering")
p+scale_color_brewer(palette="Dark2")

# Loading the Pathway Sets 
# MSigDB: http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
load(file = "../Data/reactome_genelist.RData")

# Pathway Analysis from Differential Gene Expression Data
stats <- ttop_dge$logFC
names(stats) <- ttop_dge$external_gene_name
gseaGenes <- fgseaSimple(pathways = genelist, stats = stats, nperm = 10000, 
                         minSize = 5, maxSize = Inf)

# Pathway Analysis from Differential Protein Abundance Data
stats <- ttop_prot$EGF_60_vs_PBS_60_diff
names(stats) <- ttop_prot$Gene
gseaProt <- fgseaSimple(pathways = genelist, stats = stats, nperm = 10000, 
                        minSize = 5, maxSize = Inf)

# Identifying Pathway Sets regulated on both sets (padj<=0.05)
setGenes <- gseaGenes$pathway[which(gseaGenes$padj<=0.05)]
setProt <- gseaProt$pathway[which(gseaProt$padj<=0.05)]

x <- list(setGenes = setGenes, setProt = setProt)
ggVennDiagram(x)

print(intersect(x = setGenes, y = setProt))

# Filtering DGE and DPA for common genes and retreiving p-values
data <- matrix(data = , nrow = length(common_genes), ncol = 2)
rownames(data) <- common_genes[order(common_genes)]
colnames(data) <- c("diff_genes", "diff_prot")
data[, 1] <- dge$PValue[order(dge$external_gene_name)]
data[, 2] <- prot$EGF_60_vs_PBS_60_p.val[order(prot$Gene)]
data <- as.data.frame(data)
head(data)

# Obtaining interactions from OmniPath
interactions <- import_omnipath_interactions()
interactions <- unique(as.data.frame(interactions[, 3:4]))
head(interactions)

# Transforming the obtained network into an _igraph_ object.
g <- graph_from_data_frame(d = interactions, directed = TRUE)
g <- as_graphnel(graph = g)
g

subnet <- subNetwork(rownames(data), g)
subnet

pvals <- cbind(data$diff_genes, data$diff_prot)
rownames(pvals) <- rownames(data)
pval <- aggrPvals(pvals, order = 2, plot = TRUE)

fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.5)
module <- runFastHeinz(g, scores)

logFC <- dge$logFC
names(logFC) <- dge$external_gene_name

plotModule(module, scores = scores, diff.expr = logFC)