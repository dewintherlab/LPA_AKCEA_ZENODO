#------------------------------------------------------------------------------------------
#Subset the data
ultra.str.CvD.countData <- countData[,PATIENT_TREATED]
ultra.str.CvD.rpkmData  <- rpkmData[ ,c(PATIENT_TREATED, META_COLS)]
ultra.str.CvD.colData   <- colData[   PATIENT_TREATED,]

#------------------------------------------------------------------------------------------
#Setup design and normalize the data
#Setup data
design  <- ~Response
ultra.str.CvD.dds <- DESeqDataSetFromMatrix(ultra.str.CvD.countData, ultra.str.CvD.colData, design= design)
ultra.str.CvD.dds <- DESeq(ultra.str.CvD.dds, parallel = T)

#Add extra coloring
ultra.str.CvD.colData$GR <- factor(paste(ultra.str.CvD.colData$Group,ultra.str.CvD.colData$Response, sep = "."))
colData(ultra.str.CvD.dds)$GR <- ultra.str.CvD.colData$GR
levels(ultra.str.CvD.colData$GR)

#------------------------------------------------------------------------------------------
#Call generic results, to get a clue to the variables from sample clustering

#Setup data for plotting
ultra.str.CvD.vst <- vst(ultra.str.CvD.dds)
topVarGenes <- head(order(-rowVars(assay(ultra.str.CvD.vst))),1000) #Get 1000 most variable genes
mat <- assay(ultra.str.CvD.vst)[topVarGenes,]

#Euclidean Distance plot
distsRL <- dist(t(assay(ultra.str.CvD.vst)))
d.mat <- as.matrix(distsRL)
rownames(d.mat) <- colnames(d.mat) <- paste(colnames(mat),ultra.str.CvD.colData$Flowcell, ultra.str.CvD.colData$Sex,  ultra.str.CvD.colData$Response,sep = " : ")
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("ultra.str.CvD.top1000_variable_genes.distance_plot.pdf")
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
         annotation_col = data.frame(Flowcell=ultra.str.CvD.colData$Flowcell, Sex=ultra.str.CvD.colData$Sex, Response = ultra.str.CvD.colData$Response, Patient = ultra.str.CvD.colData$Patient, row.names = row.names(ultra.str.CvD.colData)),
         annotation_colors = ann_colors,
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "ultra.str.CvD.top1000_variable_genes.kmeans_heatmap.pdf"
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(ultra.str.CvD.vst)))$x)
mat.pca$Flowcell <- ultra.str.CvD.colData$Flowcell #Add flowcell metada
mat.pca$Sex <- ultra.str.CvD.colData$Sex #Add flowcell metada
mat.pca$Patient <- ultra.str.CvD.colData$Patient #Add Patient metada
mat.pca$Response <- ultra.str.CvD.colData$Response #Add response metada
percentVar <- round(100 * prcomp(t(assay(ultra.str.CvD.vst)))$sdev^2/sum(prcomp(t(assay(ultra.str.CvD.vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, colour = Response)) +
  geom_point(size=3) +
  scale_colour_manual(values = ann_colors$Response) +
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
ggsave("ultra.str.CvD.top1000_variable_genes.PCA_plot_Response.pdf")

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
ggsave("ultra.str.CvD.top1000_variable_genes.PCA_plot_Sex.pdf")

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
ggsave("ultra.str.CvD.top1000_variable_genes.PCA_plot_Flowcell.pdf")
