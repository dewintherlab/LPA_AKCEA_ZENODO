#Call DE genes for all relevant contrasts
dir.create("AvC.differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
AvC.HEALTHY_BASELINE <- which(AvC.colData$Group == "Healthy_baseline_LPA") + 9
AvC.PATIENT_HIGH     <- which(AvC.colData$Group == "Patient_high_LPA")     + 9

RPKM.cols <- c(10:(length(colnames(AvC.countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "HEALTHY_BASELINE_v_PATIENT_HIGH"
design  <- ~Blocks + Group
dir.create(paste("AvC.differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, AvC.colData))
colnames(model.matrix(design, AvC.colData))
resultsNames(AvC.dds)

AvC.res <- results(AvC.dds, parallel= T, contrast = c("Group", "Patient_high_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("AvC.differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(AvC.res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
AvC.res <- merge(as.data.frame(AvC.res),AvC.rpkmData,by=0)

#Clean up results table
AvC.res <- AvC.res[which(!is.na(AvC.res$EntrezGene)),]
row.names(AvC.res) <- AvC.res$EntrezGene
AvC.res$EntrezGene <- NULL
colnames(AvC.res)[1] <- "RefSeq_ID"
AvC.res <- AvC.res[,c((ncol(AvC.res)-1),1,(ncol(AvC.res)),seq(from=2,to=(ncol(AvC.res)-2)))]

#Order by adjusted p value and fold change
AvC.res <- AvC.res[order(AvC.res$padj,-abs(AvC.res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(AvC.res$padj < 0.1))
length(which(AvC.res$pvalue < 0.1))

#Export the results
write.table(transform(as.data.frame(AvC.res), entrez_ID = rownames(AvC.res))[,c(length(colnames(AvC.res))+1,1:length(colnames(AvC.res)))],file=paste("AvC.differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(AvC.res, outliers = T, maxXlim = 5, minXlim = -5, autoScaleAxes = F, maxYlim = 7.5, labels = F, autoScaleLabels = T, maxLabels = 50, log2FC = 0.5)
ggsave(paste("AvC.differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
dir.create(paste("AvC.differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

#Plot the top 10 differential genes from our results, and save
for (theGene in AvC.res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(AvC.dds)), AvC.res, AvC.dds, AvC.colData)
  ggsave(paste("AvC.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#Plot the top FC genes from our results, and save
get.padj.FC.lims(AvC.res)
for (theGene in subset(AvC.res, abs(log2FoldChange) > 2 & padj < 0.1)[,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group", subgroup=1:length(colnames(AvC.dds)), AvC.res, AvC.dds, AvC.colData)
  ggsave(paste("AvC.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#------------------------------------------------------------------------------------------
#RPKM heatmap
pheatmap(subset(AvC.res, padj<0.01 & abs(log2FoldChange) > 1)[,c(AvC.HEALTHY_BASELINE, AvC.PATIENT_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(AvC.res, padj<0.01 & abs(log2FoldChange) > 1)[,"Symbol"],
         annotation_col = as.data.frame(colData(AvC.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("AvC.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.01log2FC1.pdf",  sep = "")
)

pheatmap(subset(AvC.res, padj<0.1 & abs(log2FoldChange) > 0.5)[,c(AvC.HEALTHY_BASELINE, AvC.PATIENT_HIGH)],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = as.data.frame(colData(AvC.dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 1,
         cellwidth = 20,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("AvC.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.log2FC0.5.pdf",  sep = "")
)

#Make heatmap of select genes
genes <- paste("^", c("TNF", "IL1B", "IL6", "CXCL8", "IL18", "CCL2", "CD36", "MMP2", "MMP8", "MMP9" , "MARCO" ), "$", sep = "")
AvC.res$Symbol[grep(paste(genes, collapse = "|"), AvC.res$Symbol)]

pheatmap(subset(AvC.res, Symbol %in% AvC.res$Symbol[grep(paste(genes, collapse = "|"), AvC.res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(AvC.res, Symbol %in% AvC.res$Symbol[grep(paste(genes, collapse = "|"), AvC.res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("AvC.differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = "")
)

#Check pathways
get_ontology(res = AvC.res, name = paste(contrast.name, "/GSEA", sep = ""), mydir = "AvC.differential_expression/", padj.cutoff = 0.1, log2FC.cutoff = 0.5)



