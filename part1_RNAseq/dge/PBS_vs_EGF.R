
groups <- factor(c("PBS","EGF","PBS","EGF"))


label <- c("PBS2", "EGF1", "PBS1", "EGF2")




print(colnames(subread_counts$counts))
print(colnames(subread_counts$counts))


DGEObj <- DGEList(group=groups, counts=subread_counts$counts, genes=subread_counts$annotation[,c("GeneID","Length")])
#Calculate RPKM (reads per kilobases of exon per million reads mapped) values for genes:

print(DGEObj)

design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
print(design)

contrasts <- makeContrasts(
                            PBS_vs_EGF = PBS - EGF,
                            levels=design
                          )
# restrict counts to the supplied columns
counts <- DGEObj$counts
