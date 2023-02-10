#Call DE genes for all relevant contrasts
dir.create("AvB.differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
AvB.HEALTHY_BASELINE <- which(AvB.colData$Group == "Healthy_baseline_LPA") + 9
AvB.HEALTHY_HIGH     <- which(AvB.colData$Group == "Healthy_high_LPA")     + 9

RPKM.cols <- c(10:(length(colnames(AvB.countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "HEALTHY_BASELINE_v_HIGH"
design  <- ~Group
dir.create(paste("AvB.differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, AvB.colData))
colnames(model.matrix(design, AvB.colData))
resultsNames(AvB.dds)

AvB.res <- results(AvB.dds, parallel= T, contrast = c("Group", "Healthy_high_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("AvB.differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(AvB.res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
AvB.res <- merge(as.data.frame(AvB.res),AvB.rpkmData,by=0)

#Clean up results table
AvB.res <- AvB.res[which(!is.na(AvB.res$EntrezGene)),]
row.names(AvB.res) <- AvB.res$EntrezGene
AvB.res$EntrezGene <- NULL
colnames(AvB.res)[1] <- "RefSeq_ID"
AvB.res <- AvB.res[,c((ncol(AvB.res)-1),1,(ncol(AvB.res)),seq(from=2,to=(ncol(AvB.res)-2)))]

#Order by adjusted p value and fold change
AvB.res <- AvB.res[order(AvB.res$padj,-abs(AvB.res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(AvB.res$padj < 0.1))
length(which(AvB.res$pvalue < 0.1))

#Export the results
write.table(transform(as.data.frame(AvB.res), entrez_ID = rownames(AvB.res))[,c(length(colnames(AvB.res))+1,1:length(colnames(AvB.res)))],file=paste("AvB.differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(AvB.res, outliers = T, maxXlim = 3, minXlim = -3, autoScaleAxes = F, maxYlim = 3, labels = T, autoScaleLabels = T, maxLabels = 10, padj = 0.1, log2FC = 0.5)
ggsave(paste("AvB.differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot_more_zoom_no_labels.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
dir.create(paste("AvB.differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

#Plot the top 10 differential genes from our results, and save
for (theGene in AvB.res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(AvB.dds)), AvB.res, AvB.dds, AvB.colData)
  ggsave(paste("AvB.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#Plot the top FC genes from our results, and save
for (theGene in subset(AvB.res, abs(log2FoldChange) > 0.5 & padj < 0.1)[,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(AvB.dds)), AvB.res, AvB.dds, AvB.colData)
  ggsave(paste("AvB.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#------------------------------------------------------------------------------------------
#RPKM heatmap
pheatmap(subset(AvB.res, padj<0.1 & abs(log2FoldChange) > 0.5)[,c(AvB.HEALTHY_BASELINE, AvB.HEALTHY_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(AvB.res, padj<0.1 & abs(log2FoldChange) > 0.5)[,"Symbol"],
         annotation_col = as.data.frame(colData(AvB.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("AvB.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1log2FC0.5.pdf",  sep = "")
)

#Make heatmap of select genes
genes <- paste("^", c("TNF", "IL1B", "IL6", "CXCL8", "IL18", "CCL2", "CD36", "MMP2", "MMP8", "MMP9" , "MARCO" ), "$", sep = "")
AvB.res$Symbol[grep(paste(genes, collapse = "|"), AvB.res$Symbol)]

pheatmap(subset(AvB.res, Symbol %in% AvB.res$Symbol[grep(paste(genes, collapse = "|"), AvB.res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(AvB.res, Symbol %in% AvB.res$Symbol[grep(paste(genes, collapse = "|"), AvB.res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("AvB.differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = "")
         )

#Check pathways
get_ontology(res = AvB.res, name = paste(contrast.name, "/GSEA", sep = ""), mydir = "AvB.differential_expression/", padj.cutoff = 0.25)


#replot some GSEA ES plots
replotGSEA(path = "/Users/koen/gsea_home/output/may02/AvC_REACTOME_IFN_GAMMA_SIGNALLING.Gsea.1556804607719/", gene.set = "REACTOME", class.name = "Healthy Control")

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#PCA plot of A B, and C
ABC.countData <- countData[,c(HEALTHY_BASELINE, HEALTHY_HIGH, PATIENT_HIGH)]
ABC.colData   <- colData[   c(HEALTHY_BASELINE, HEALTHY_HIGH, PATIENT_HIGH),]

#------------------------------------------------------------------------------------------
#Setup design and normalize the data
#Setup data
design  <- ~Group
ABC.dds <- DESeqDataSetFromMatrix(ABC.countData, ABC.colData, design= design)
ABC.dds <- DESeq(ABC.dds, parallel = T)

#------------------------------------------------------------------------------------------
#Call generic results, to get a clue to the variables from sample clustering

#Setup data for plotting
ABC.vst <- vst(ABC.dds)
IFN.genes <- scan(file="REACTOME_INTERFERON_GAMMA_SIGNALING.grp", what = "list")
IFN.genes <- IFN.genes[9:length(IFN.genes)]
IFN.genes <- rownames(rpkmData[rpkmData$Symbol %in% IFN.genes,])
topVarGenes <- which(rownames(assay(ABC.vst)) %in% IFN.genes) #Get All IFN pathway genes

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(ABC.vst)[topVarGenes,]))$x)
mat.pca$Group <- ABC.colData$Group #Add group metadata
mat.pca$Flowcell <- ABC.colData$Flowcell #Add flowcell metada
mat.pca$Sex <- ABC.colData$Sex #Add flowcell metada
percentVar <- round(100 * prcomp(t(assay(ABC.vst)[topVarGenes,]))$sdev^2/sum(prcomp(t(assay(ABC.vst)[topVarGenes,]))$sdev^2))

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
        legend.position = "right",
        aspect.ratio = 1:1
  )
ggsave("ABC.IFN_genes.PCA_plot_Group.pdf")

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
        legend.position = "right",
        aspect.ratio = 1:1
  )
ggsave("ABC.top1000_variable_genes.PCA_plot_Sex.pdf")

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
        legend.position = "right",
        aspect.ratio = 1:1
  )
ggsave("ABC.top1000_variable_genes.PCA_plot_Flowcell.pdf")


#PRKM box ABC
ATH.genes <- scan(file="ATHEROGENIC_GENES.txt", what = "list")
ATH.genes <- ATH.genes[9:length(ATH.genes)]
AvC.sig.res <- subset(AvC.res, padj < 0.1 & abs(log2FoldChange) > 0.5)

topVarGenes <- rpkmData$Symbol %in% AvC.sig.res[AvC.sig.res$Symbol %in% ATH.genes,"Symbol"]
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Healthy Baseline","Healthy High LPA")
topdiff.res <- as.matrix(topdiff.res)

topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility
m <- melt(topdiff.res)
m <- melt(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Healthy Baseline"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m$Symbol <- rpkmData[topVarGenes,][rownames(rpkmData[topVarGenes,]) %in% m$refSeq_ID,"Symbol"]
m <- m[m$RPKM < 10,]

#Custom jitter function - call before each ggplot call where needed
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

ggplot(m,aes(x = Group, y = RPKM, col = Group)) +
  geom_point(position=myjit,size=3,alpha = 0.25) +
  #scale_y_log10() +
  ylab("RPKM") +
  scale_x_discrete(limits = colnames(topdiff.res)[c(2,3,1)]) +
  scale_colour_manual(name="Group", values = c("firebrick1","chartreuse3","coral"), breaks = colnames(topdiff.res), labels = c("Healthy controls - Baseline LPA level",
                                                                                                             "Healthy controls - High LPA level",
                                                                                                             "Patients - High LPA level")) +
  ggtitle("Median RPKM of differentially expressed atherogenic genes") +
  stat_summary(fun.y = median,
               fun.ymin = function(z) {quantile(z, 0.25)},
               fun.ymax = function(z) {quantile(z, 0.75)},
               geom = "crossbar",
               width = 0.5,
               size = 1,
               mapping = aes( color = Group )
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z,0.25)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
               fun.ymax = function(z) {quantile(z,0.95)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z, 0.05)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
               fun.ymax = function(z) {quantile(z, 0.95)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  geom_label_repel(
    aes(label = as.character(Symbol)),
    box.padding = 0.35, 
    position = myjit,
    point.padding = 0.5
   ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "none")
ggsave("ABC.ATHgenes.Median_RPKM_increased_plot.pdf")


#PRKM trend ABCD
ATH.genes <- scan(file="ATHEROGENIC_GENES.txt", what = "list")
ATH.genes <- ATH.genes[9:length(ATH.genes)]
AvC.sig.res <- subset(AvC.res, padj < 0.1 & abs(log2FoldChange) > 0.5)

topVarGenes <- rpkmData$Symbol %in% AvC.sig.res[AvC.sig.res$Symbol %in% ATH.genes,"Symbol"]
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Healthy Baseline","Healthy High LPA","Patient Treated")
topdiff.res <- as.matrix(topdiff.res)

topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility
m <- melt(topdiff.res)
m <- melt(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Healthy Baseline"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m$Symbol <- rpkmData[topVarGenes,][rownames(rpkmData[topVarGenes,]) %in% m$refSeq_ID,"Symbol"]
m <- m[m$RPKM < 5,]

#Custom jitter function - call before each ggplot call where needed
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

ggplot(m,aes(x = Group, y = RPKM, col = Group)) +
  #scale_y_log10() +
  ylab("RPKM") +
  scale_x_discrete(limits = colnames(topdiff.res)[c(2,3,1,4)],expand = c(0,0)) +
  scale_y_continuous(limits= c(0,5),expand = c(0,0)) +
  scale_colour_manual(name="Group", values = c("firebrick1","chartreuse3","coral","chartreuse4"), breaks = colnames(topdiff.res), labels = c("Healthy controls - Baseline LPA level",
                                                                                                                               "Healthy controls - High LPA level",
                                                                                                                               "Patients - High LPA level", "Patients - Treated")) +
  ggtitle("Median RPKM of differentially expressed atherogenic genes") +
  geom_smooth(se = T, 
              method = "loess", 
              stat = "summary", 
              fun.y = median,
              fun.ymin = function(z) {quantile(z, 0.25)},
              fun.ymax = function(z) {quantile(z, 0.75)},  
              mapping = aes(group = 1, y = RPKM), 
              color = "dodgerblue",
              fill = "lightblue1",
              size = 3
  ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "none")
ggsave("ABCD.ATHgenes.Median_RPKM_increased_plot.pdf")

#PRKM box ABCD
ATH.genes <- scan(file="ATHEROGENIC_GENES.txt", what = "list")
ATH.genes <- ATH.genes[9:length(ATH.genes)]
AvC.sig.res <- subset(AvC.res, padj < 0.1 & abs(log2FoldChange) > 0.5)

topVarGenes <- rpkmData$Symbol %in% AvC.sig.res[AvC.sig.res$Symbol %in% ATH.genes,"Symbol"]
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Healthy Baseline","Healthy High LPA","Patient Treated")
topdiff.res <- as.matrix(topdiff.res)

topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility
m <- melt(topdiff.res)
m <- melt(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Healthy Baseline"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m$Symbol <- rpkmData[topVarGenes,][rownames(rpkmData[topVarGenes,]) %in% m$refSeq_ID,"Symbol"]
m <- m[m$RPKM < 5,]

#Custom jitter function - call before each ggplot call where needed
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

ggplot(m,aes(x = Group, y = RPKM, col = Group)) +
  #scale_y_log10() +
  ylab("RPKM") +
  scale_x_discrete(limits = colnames(topdiff.res)[c(2,3,1,4)]) +
  scale_y_continuous(limits= c(0,5)) +
  scale_colour_manual(name="Group", values = c("firebrick1","chartreuse3","coral","chartreuse4"), breaks = colnames(topdiff.res), labels = c("Healthy controls - Baseline LPA level",
                                                                                                                                             "Healthy controls - High LPA level",
                                                                                                                                             "Patients - High LPA level", "Patients - Treated")) +
  ggtitle("Median RPKM of differentially expressed atherogenic genes") +
  geom_point(position=myjit,size=3,alpha = 0.25) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) {quantile(z, 0.25)},
               fun.ymax = function(z) {quantile(z, 0.75)},
               geom = "crossbar",
               width = 0.5,
               size = 1,
               mapping = aes( color = Group )
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z,0.25)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
               fun.ymax = function(z) {quantile(z,0.95)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z, 0.05)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
               fun.ymax = function(z) {quantile(z, 0.95)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  geom_label_repel(
    aes(label = as.character(Symbol)),
    box.padding = 0.35, 
    position = myjit,
    point.padding = 0.5
  ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "none")
ggsave("ABCD.ATHgenes.Median_RPKM_increased_boxplot.pdf")



#---------------------------------------------------------------------
#Adding results from ANITSCHKOW (aPCSK9) study
#Read tables
aPCSK9.rpkmData <- read.table(file = "../LPA_PCSK9/RPKM.table",     header = T, row.names = 1)
aPCSK9.colData  <- read.table(file = "../LPA_PCSK9/metadata.table", header = T, row.names = 1, colClasses = "factor")

#Define columns
PCSK9_BASELINE   <- which(aPCSK9.colData$Group == "BASELINE")
PCSK9_EVOLOCUMAB <- which(aPCSK9.colData$Group == "EVOLOCUMAB")
PCSK9_PLACEBO    <- which(aPCSK9.colData$Group == "PLACEBO")

#Subset genes we want to plot from 'common' RPKM table
ATH.genes <- scan(file="ATHEROGENIC_GENES.txt", what = "list")
ATH.genes <- ATH.genes[9:length(ATH.genes)]
AvC.sig.res <- subset(AvC.res, padj < 0.1 & abs(log2FoldChange) > 0.5)
topVarGenes <- rpkmData$Symbol %in% AvC.sig.res[AvC.sig.res$Symbol %in% ATH.genes,"Symbol"]
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED],    1,median))
colnames(topdiff.res) <- c("Patient High LPA - AKCEA","Healthy Baseline","Healthy High LPA","Patient Treated - LP(a) siRNA")
topdiff.res <- as.matrix(topdiff.res)
topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility

#Add ANITSCHKOW data
aPCSK9_subset.rpkmData <- aPCSK9.rpkmData[row.names(aPCSK9.rpkmData) %in% row.names(topdiff.res),]
aPCSK9_subset.res      <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_BASELINE],1,median))
aPCSK9_subset.res[,2]  <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_EVOLOCUMAB],1,median))
aPCSK9_subset.res[,3]  <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_PLACEBO],1,median))
colnames(aPCSK9_subset.res) <- c("Patient High LPA - ANITSCHKOW","Patient Treated - anti-PCSK9 Ab.", "Patient Treated - Placebo")
aPCSK9_subset.res <- as.matrix(aPCSK9_subset.res)

merged.topdiff.res <- merge(topdiff.res, aPCSK9_subset.res, by = 0)
row.names(merged.topdiff.res) <- merged.topdiff.res$Row.names
merged.topdiff.res$Row.names  <- NULL
merged.topdiff.res <- as.matrix(merged.topdiff.res)


#Prep for plot
m <- melt(merged.topdiff.res)
m <- melt(merged.topdiff.res[which(merged.topdiff.res[,"Patient High LPA - AKCEA"] > merged.topdiff.res[,"Healthy Baseline"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m <- merge(m, rpkmData[row.names(rpkmData) %in% m$refSeq_ID,"Symbol", drop = F], by.x = 1, by.y = 0)
m <- m[m$RPKM < 5,]

#Custom jitter function - call before each ggplot call where needed
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

ggplot(m,aes(x = Group, y = RPKM, col = Group)) +
  #scale_y_log10() +
  ylab("RPKM") +
  scale_x_discrete(limits = colnames(merged.topdiff.res)[c(2,3,1,5,4,6,7)]) +
  scale_y_continuous(limits= c(0,5)) +
  scale_colour_manual(name="Group", values = c("firebrick1","chartreuse3","coral","chartreuse4","darkblue","blue","dodgerblue"), breaks = colnames(merged.topdiff.res), labels = c("Healthy controls - Baseline LPA level",
                                                                                                                                             "Healthy controls - High LPA level",
                                                                                                                                             "Patients - High LPA level (AKCEA)",
                                                                                                                                             "Patients - High LPA level (ANITSCHKOW)",
                                                                                                                                             "Patients - Treated (LP(a) siRNA)",
                                                                                                                                             "Patients - Treated (anti-PCSK9 Ab.)",
                                                                                                                                             "Patients - Treated (Placebo)")) +
  ggtitle("Median RPKM of differentially expressed atherogenic genes") +
  geom_point(position=myjit,size=3,alpha = 0.25) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) {quantile(z, 0.25)},
               fun.ymax = function(z) {quantile(z, 0.75)},
               geom = "crossbar",
               width = 0.5,
               size = 1,
               mapping = aes( color = Group )
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z,0.25)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
               fun.ymax = function(z) {quantile(z,0.95)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z, 0.05)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
               fun.ymax = function(z) {quantile(z, 0.95)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "none", 
        aspect.ratio = 2.5/3)
ggsave("ABCDANIT.ATHgenes.Median_RPKM_increased_boxplot.pdf")



#--------------------------------------------------------------------------
#Let's compare all AvC genes

#Subset genes we want to plot from 'common' RPKM table
AvC.sig.res <- subset(AvC.res, padj < 0.1 & abs(log2FoldChange) > 0.5)
topVarGenes <- row.names(rpkmData[rpkmData$Symbol %in% AvC.sig.res$Symbol,])
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED],    1,median))
colnames(topdiff.res) <- c("Patient High LPA - AKCEA","Healthy Baseline","Healthy High LPA","Patient Treated - LP(a) siRNA")
topdiff.res <- as.matrix(topdiff.res)
#topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility

#Add ANITSCHKOW data
aPCSK9_subset.rpkmData <- aPCSK9.rpkmData[row.names(aPCSK9.rpkmData) %in% row.names(topdiff.res),]
aPCSK9_subset.res      <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_BASELINE],1,median))
aPCSK9_subset.res[,2]  <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_EVOLOCUMAB],1,median))
aPCSK9_subset.res[,3]  <- as.data.frame(apply(aPCSK9_subset.rpkmData[,PCSK9_PLACEBO],1,median))
colnames(aPCSK9_subset.res) <- c("Patient High LPA - ANITSCHKOW","Patient Treated - anti-PCSK9 Ab.", "Patient Treated - Placebo")
aPCSK9_subset.res <- as.matrix(aPCSK9_subset.res)

merged.topdiff.res <- merge(topdiff.res, aPCSK9_subset.res, by = 0)
row.names(merged.topdiff.res) <- merged.topdiff.res$Row.names
merged.topdiff.res$Row.names  <- NULL
merged.topdiff.res <- as.matrix(merged.topdiff.res)


#Prep for plot
m <- melt(merged.topdiff.res)
m <- melt(merged.topdiff.res[which(merged.topdiff.res[,"Patient High LPA - AKCEA"] > merged.topdiff.res[,"Healthy Baseline"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m <- merge(m, rpkmData[row.names(rpkmData) %in% m$refSeq_ID,"Symbol", drop = F], by.x = 1, by.y = 0)
m <- m[m$RPKM < 250,]

#Custom jitter function - call before each ggplot call where needed
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.3,
                 dodge.width = 0.1,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )

ggplot(m,aes(x = Group, y = RPKM, col = Group)) +
  #scale_y_log10() +
  ylab("RPKM") +
  scale_x_discrete(limits = colnames(merged.topdiff.res)[c(2,3,1,5,4,6,7)]) +
  #scale_y_continuous(limits= c(0,5)) +
  scale_colour_manual(name="Group", values = c("firebrick1","chartreuse3","coral","chartreuse4","darkblue","blue","dodgerblue"), breaks = colnames(merged.topdiff.res), labels = c("Healthy controls - Baseline LPA level",
                                                                                                                                                                                   "Healthy controls - High LPA level",
                                                                                                                                                                                   "Patients - High LPA level (AKCEA)",
                                                                                                                                                                                   "Patients - High LPA level (ANITSCHKOW)",
                                                                                                                                                                                   "Patients - Treated (LP(a) siRNA)",
                                                                                                                                                                                   "Patients - Treated (anti-PCSK9 Ab.)",
                                                                                                                                                                                   "Patients - Treated (Placebo)")) +
  ggtitle("Median RPKM of differentially expressed AvC genes higher in C") +
  geom_point(position=myjit,size=3,alpha = 0.25) +
  stat_summary(fun.y = median,
               fun.ymin = function(z) {quantile(z, 0.25)},
               fun.ymax = function(z) {quantile(z, 0.75)},
               geom = "crossbar",
               width = 0.5,
               size = 1,
               mapping = aes( color = Group )
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z,0.25)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
               fun.ymax = function(z) {quantile(z,0.95)},
               geom = "linerange",
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
               fun.ymax = function(z) {quantile(z, 0.05)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
               fun.ymax = function(z) {quantile(z, 0.95)},
               geom = "errorbar",
               width = 0.25,
               size = 0.5,
               mapping = aes(color = Group)
  ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "none", 
        aspect.ratio = 2.5/3)
ggsave("ABCDANIT.AvC_higher_in_C_under_250.Median_RPKM_increased_boxplot.pdf")
