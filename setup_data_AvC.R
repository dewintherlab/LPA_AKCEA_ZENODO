#------------------------------------------------------------------------------------------
#Subset the data
AvC.countData <- countData[,c(HEALTHY_BASELINE, PATIENT_HIGH)]
AvC.rpkmData  <- rpkmData[ ,c(HEALTHY_BASELINE, PATIENT_HIGH, META_COLS)]
AvC.colData   <- colData[   c(HEALTHY_BASELINE, PATIENT_HIGH),]

#------------------------------------------------------------------------------------------
#Setup design and normalize the data
#Setup data
design  <- ~Blocks + Group
AvC.dds <- DESeqDataSetFromMatrix(AvC.countData, AvC.colData, design= design)
AvC.dds <- DESeq(AvC.dds, parallel = T)

#------------------------------------------------------------------------------------------
#Call generic results, to get a clue to the variables from sample clustering

#Setup data for plotting
AvC.vst <- vst(AvC.dds)
topVarGenes <- head(order(-rowVars(assay(AvC.vst))),1000) #Get 1000 most variable genes
mat <- assay(AvC.vst)[topVarGenes,]

#Euclidean Distance plot
distsRL <- dist(t(assay(AvC.vst)))
d.mat <- as.matrix(distsRL)
rownames(d.mat) <- colnames(d.mat) <- paste(colnames(mat),AvC.colData$Flowcell, AvC.colData$Sex, AvC.colData$Group,sep = " : ")
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("AvC.top1000_variable_genes.distance_plot.pdf")
heatmap.2(d.mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13),cexRow = 0.5,cexCol = 0.5)
dev.off()

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = data.frame(Flowcell=AvC.colData$Flowcell, Sex=AvC.colData$Sex, Group=AvC.colData$Group, row.names = row.names(AvC.colData)),
         annotation_colors = ann_colors,
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "AvC.top1000_variable_genes.kmeans_heatmap.pdf"
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(AvC.vst)))$x)
mat.pca$Group <- AvC.colData$Group #Add group metadata
mat.pca$Flowcell <- AvC.colData$Flowcell #Add flowcell metada
mat.pca$Sex <- AvC.colData$Sex #Add flowcell metada
percentVar <- round(100 * prcomp(t(assay(AvC.vst)))$sdev^2/sum(prcomp(t(assay(AvC.vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Group)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Group) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("AvC.top1000_variable_genes.PCA_plot_Group.pdf")

ggplot(mat.pca, aes(PC1, PC2, color = Sex)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Sex) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("AvC.top1000_variable_genes.PCA_plot_Sex.pdf")

ggplot(mat.pca, aes(PC1, PC2, color = Flowcell)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar[2],"% variance")) +
  ylab(paste0("PC2: ",percentVar[3],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("AvC.top1000_variable_genes.PCA_plot_Flowcell.pdf")

