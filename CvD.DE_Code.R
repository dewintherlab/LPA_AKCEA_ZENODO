#Call DE genes for all relevant contrasts
dir.create("CvD.differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
CvD.PATIENT_TREATED <- which(CvD.colData$Group == "Patient_treated_LPA") + 9
CvD.PATIENT_HIGH     <- which(CvD.colData$Group == "Patient_high_LPA")     + 9

RPKM.cols <- c(10:(length(colnames(CvD.countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "PATIENT_TREATED_v_PATIENT_HIGH"
design  <- ~Patient + Group
dir.create(paste("CvD.differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, CvD.colData))
colnames(model.matrix(design, CvD.colData))
resultsNames(CvD.dds)

CvD.res <- results(CvD.dds, parallel= T, contrast = c("Group", "Patient_treated_LPA", "Patient_high_LPA"))

#Visualise
pdf(paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(CvD.res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
CvD.res <- merge(as.data.frame(CvD.res),CvD.rpkmData,by=0)

#Clean up results table
CvD.res <- CvD.res[which(!is.na(CvD.res$EntrezGene)),]
row.names(CvD.res) <- CvD.res$EntrezGene
CvD.res$EntrezGene <- NULL
colnames(CvD.res)[1] <- "RefSeq_ID"
CvD.res <- CvD.res[,c((ncol(CvD.res)-1),1,(ncol(CvD.res)),seq(from=2,to=(ncol(CvD.res)-2)))]

#Order by adjusted p value and fold change
CvD.res <- CvD.res[order(CvD.res$padj,-abs(CvD.res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(CvD.res$padj < 0.1))
length(which(CvD.res$pvalue < 0.1))

#Export the results
write.table(transform(as.data.frame(CvD.res), entrez_ID = rownames(CvD.res))[,c(length(colnames(CvD.res))+1,1:length(colnames(CvD.res)))],file=paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(CvD.res, outliers = T, maxXlim = 2, minXlim = -2, autoScaleAxes = F, maxYlim = 7.5, labels = F, autoScaleLabels = F, maxLabels = 50, log2FC = 0.5)
ggsave(paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
dir.create(paste("CvD.differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

#Plot the top 10 differential genes from our results, and save
for (theGene in CvD.res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(CvD.dds)), CvD.res, CvD.dds, CvD.colData)
  ggsave(paste("CvD.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#Plot the top FC genes from our results, and save
get.padj.FC.lims(CvD.res)
for (theGene in subset(CvD.res, abs(log2FoldChange) > 0.4 & padj < 0.1)[,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(CvD.dds)), CvD.res, CvD.dds, CvD.colData)
  ggsave(paste("CvD.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#Plot some interesting genes
for (theGene in c("CCR2", "CX3CR1", "TLR2", "TLR4", "MYD88")){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(CvD.dds)), CvD.res, CvD.dds, CvD.colData)
  ggsave(paste("CvD.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#------------------------------------------------------------------------------------------
#RPKM heatmap
pheatmap(subset(CvD.res, padj<0.1 & abs(log2FoldChange) > 0.5)[,c(CvD.PATIENT_TREATED, CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(CvD.res, padj<0.1 & abs(log2FoldChange) > 0.5)[,"Symbol"],
         annotation_col = as.data.frame(colData(CvD.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.log2FC0.5.pdf",  sep = "")
)

pheatmap(subset(CvD.res, padj<0.1 & abs(log2FoldChange) > 0.1)[,c(CvD.PATIENT_TREATED, CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = as.data.frame(colData(CvD.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 1,
         cellwidth = 20,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.log2FC0.1.pdf",  sep = "")
)

pheatmap(subset(CvD.res, padj<0.1)[,c(CvD.PATIENT_TREATED, CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = as.data.frame(colData(CvD.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 1,
         cellwidth = 20,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.pdf",  sep = "")
)
#Make heatmap of select genes
genes <- paste("^", c("TNF", "IL1B", "IL6", "CXCL8", "IL18", "CCL2", "CD36", "MMP2", "MMP8", "MMP9" , "MARCO" ), "$", sep = "")
CvD.res$Symbol[grep(paste(genes, collapse = "|"), CvD.res$Symbol)]

pheatmap(subset(CvD.res, Symbol %in% CvD.res$Symbol[grep(paste(genes, collapse = "|"), CvD.res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(CvD.res, Symbol %in% CvD.res$Symbol[grep(paste(genes, collapse = "|"), CvD.res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("CvD.differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = "")
)
#And plot them
dir.create(paste("CvD.differential_expression/", contrast.name, "/select INF gene counts/", sep = ""), showWarnings = FALSE)
listofgenes <- CvD.res$Symbol[grep(paste(genes, collapse = "|"), CvD.res$Symbol)]
for (theGene in listofgenes){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(CvD.dds)), CvD.res, CvD.dds, CvD.colData)
  ggsave(paste("CvD.differential_expression/", contrast.name, "/select INF gene counts/", theGene, " counts.pdf", sep = ""))
}


#Check pathways
get_ontology(res = CvD.res, name = paste(contrast.name, "/GSEA", sep = ""), mydir = "CvD.differential_expression/", padj.cutoff = 0.1)



