library("DEP")
library("readr")
library("vsn")
library("dplyr")
library("tidyr")
library("limma")
library("ggplot2")
library("ggrepel")
library("knitr")
library("tidyverse")

data <- read.delim("../Data/protein_intensities.txt")
colnames(data)
head(data)

# Are there any duplicated protein names?
data$Accession %>% duplicated() %>% any()

# Make a table of duplicated protein names
data %>% group_by(Accession) %>% summarize(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Accession" and "Sequence" column.
data_unique <- make_unique(data, "Accession", "Sequence", delim = ";")
head(data_unique)

# Are there any duplicated names?
data$name %>% duplicated() %>% any()

# We show the sample ID's
labels <- colnames(data_unique)[4:43]
print(labels)

# We setup the experimental design matrix
experimental_design <- matrix(data = , nrow = length(labels), ncol = 3)
experimental_design[, 1] <- labels
experimental_design[, 2] <- sapply(strsplit(x = labels, split = "_Rep", fixed = TRUE), "[", 1)
experimental_design[, 3] <- sapply(strsplit(x = labels, split = "_", fixed = TRUE), "[", 3)
colnames(experimental_design) <- c("label", "condition", "replicate")
experimental_design <- as.data.frame(experimental_design)
print(experimental_design)

# Generate a SummarizedExperiment object using an experimental design
Int_columns <- 4:43 # get Intensity column numbers
data_se <- make_se(data_unique, Int_columns, as.data.frame(experimental_design))
data_se

# Filter for proteins that are identified in 2 out of 4 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 2)
data_filt

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Percentage of missing data
perc <- length(which(is.na(data_norm@assays@data@listData[[1]])))/(dim(data_norm@assays@data@listData[[1]])[1]*dim(data_norm@assays@data@listData[[1]])[2])
print(paste0(round(x = 100*perc, digits = 2), "%"))

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# We can then visualize the effect of imputation on the data
plot_imputation(data_norm, data_imp)

matrix_data <- data_imp@assays@data@listData[[1]]

# PCA plot based on condition
groups <- experimental_design$condition
names(groups) <- experimental_design$label

data.pca <- t(matrix_data)
data.pca <- cbind(data.pca, as.matrix(groups))
colnames(data.pca)[ncol(data.pca)] <- "Group"
data.pca <- as.data.frame(data.pca)
data.pca[, 1:(ncol(data.pca)-1)] <- lapply(data.pca[, 1:(ncol(data.pca)-1)], 
                                           function(x) as.numeric(as.character(x)))
res.pca <- prcomp(data.pca[, -ncol(data.pca)], scale. = TRUE)
res.plot <-  as.data.frame(cbind(res.pca$x[, 1], res.pca$x[, 2], 
                                 as.character(data.pca$Group), rownames(data.pca)))
res.plot[, 1:2] <- lapply(res.plot[, 1:2], function(x) as.numeric(as.character(x)))
res.plot[, 3:4] <- lapply(res.plot[, 3:4], function(x) as.character(x))
colnames(res.plot) <- c("pc1", "pc2", "Group", "sample")
percentages <- ((res.pca$sdev)^2 / sum(res.pca$sdev^2)*100)[1:2]
pp <- ggplot(res.plot, aes(x=pc1, y=pc2, color=Group)) +
  geom_point(size=7, alpha = 0.5) +
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(x = percentages[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(x = percentages[2], digits = 2), "%)")) +
  xlim(c(-max(abs(res.pca$x[, 1])),max(abs(res.pca$x[, 1])))) +
  ylim(c(-max(abs(res.pca$x[, 2])),max(abs(res.pca$x[, 2])))) + 
  theme(legend.position = "none") +
  geom_text_repel(data = res.plot, aes(label=sample), max.overlaps = 100)
plot(pp)

# PCA plot based on relpicate
groups <- experimental_design$replicate
names(groups) <- experimental_design$label

data.pca <- t(matrix_data)
data.pca <- cbind(data.pca, as.matrix(groups))
colnames(data.pca)[ncol(data.pca)] <- "Group"
data.pca <- as.data.frame(data.pca)
data.pca[, 1:(ncol(data.pca)-1)] <- lapply(data.pca[, 1:(ncol(data.pca)-1)], 
                                           function(x) as.numeric(as.character(x)))
res.pca <- prcomp(data.pca[, -ncol(data.pca)], scale. = TRUE)
res.plot <-  as.data.frame(cbind(res.pca$x[, 1], res.pca$x[, 2], 
                                 as.character(data.pca$Group), rownames(data.pca)))
res.plot[, 1:2] <- lapply(res.plot[, 1:2], function(x) as.numeric(as.character(x)))
res.plot[, 3:4] <- lapply(res.plot[, 3:4], function(x) as.character(x))
colnames(res.plot) <- c("pc1", "pc2", "Group", "sample")
percentages <- ((res.pca$sdev)^2 / sum(res.pca$sdev^2)*100)[1:2]
pp <- ggplot(res.plot, aes(x=pc1, y=pc2, color=Group)) +
  geom_point(size=7, alpha = 0.5) +
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(x = percentages[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(x = percentages[2], digits = 2), "%)")) +
  xlim(c(-max(abs(res.pca$x[, 1])),max(abs(res.pca$x[, 1])))) +
  ylim(c(-max(abs(res.pca$x[, 2])),max(abs(res.pca$x[, 2])))) + 
  theme(legend.position = "none") +
  geom_text_repel(data = res.plot, aes(label=sample), max.overlaps = 100)
plot(pp)

# Remove Batch effects
batch_rem <- removeBatchEffect(x = matrix_data, batch = experimental_design$replicate)

# PCA plot based on condition after removing the batch effects
groups <- experimental_design$condition
names(groups) <- experimental_design$label

data.pca <- t(batch_rem)
data.pca <- cbind(data.pca, as.matrix(groups))
colnames(data.pca)[ncol(data.pca)] <- "Group"
data.pca <- as.data.frame(data.pca)
data.pca[, 1:(ncol(data.pca)-1)] <- lapply(data.pca[, 1:(ncol(data.pca)-1)], 
                                           function(x) as.numeric(as.character(x)))
res.pca <- prcomp(data.pca[, -ncol(data.pca)], scale. = TRUE)
res.plot <-  as.data.frame(cbind(res.pca$x[, 1], res.pca$x[, 2], 
                                 as.character(data.pca$Group), rownames(data.pca)))
res.plot[, 1:2] <- lapply(res.plot[, 1:2], function(x) as.numeric(as.character(x)))
res.plot[, 3:4] <- lapply(res.plot[, 3:4], function(x) as.character(x))
colnames(res.plot) <- c("pc1", "pc2", "Group", "sample")
percentages <- ((res.pca$sdev)^2 / sum(res.pca$sdev^2)*100)[1:2]
pp <- ggplot(res.plot, aes(x=pc1, y=pc2, color=Group)) +
  geom_point(size=7, alpha = 0.5) +
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  #geom_path(arrow=arrow()) +
  theme_minimal() +
  xlab(paste0("PC1 (", round(x = percentages[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(x = percentages[2], digits = 2), "%)")) +
  xlim(c(-max(abs(res.pca$x[, 1])),max(abs(res.pca$x[, 1])))) +
  ylim(c(-max(abs(res.pca$x[, 2])),max(abs(res.pca$x[, 2])))) + 
  theme(legend.position = "none") +
  geom_text_repel(data = res.plot, aes(label=sample), max.overlaps = 100)
plot(pp)

# Create the DEP object after removing the batch effects
data_batch <- data_imp
data_batch@assays@data@listData[[1]] <- batch_rem

# Plot intensity distributions before and after batch correction
plot_normalization(data_imp, data_batch)

# Test manually defined comparisons for time-point 60.
data_diff_manual <- test_diff(data_batch, type = "manual", 
                              test = c("EGF_60_vs_PBS_60"))

ttop <- as.data.frame(data_diff_manual@elementMetadata)
ttop <- ttop[, c(1:4, 10:12)]
head(ttop)

ttop$expression = ifelse(ttop$EGF_60_vs_PBS_60_p.val < 0.05 & abs(ttop$EGF_60_vs_PBS_60_diff) >= 1, 
                         ifelse(ttop$EGF_60_vs_PBS_60_diff> 1 ,'Up','Down'),
                         'Stable')

ttop$label <- NA
ttop$label[which(ttop$expression!="Stable")] <- ttop$Gene[which(ttop$expression!="Stable")]

p <- ggplot(data=ttop, aes(x=EGF_60_vs_PBS_60_diff, y=-log10(EGF_60_vs_PBS_60_p.val), col=expression, label=label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential Protein Abundance")
p