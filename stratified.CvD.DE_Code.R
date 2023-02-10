#Call DE genes for all relevant contrasts
dir.create("str.CvD.differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
str.CvD.PATIENT_TREATED <- which(str.CvD.colData$Group == "Patient_treated_LPA") + 9
str.CvD.PATIENT_HIGH     <- which(str.CvD.colData$Group == "Patient_high_LPA")     + 9

RPKM.cols <- c(10:(length(colnames(str.CvD.countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "PATIENT_HIGH_RESPONDER_v_PATIENT_LOW_RESPONDER"
design  <- ~Group * Response
dir.create(paste("str.CvD.differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, str.CvD.colData))
colnames(model.matrix(design, str.CvD.colData))
resultsNames(str.CvD.dds)

str.CvD.res <- results(str.CvD.dds, parallel= T, contrast = list(c("Response_High_vs_Low", "GroupPatient_treated_LPA.ResponseHigh")))

#Visualise
pdf(paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(str.CvD.res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
str.CvD.res <- merge(as.data.frame(str.CvD.res),str.CvD.rpkmData,by=0)

#Clean up results table
str.CvD.res <- str.CvD.res[which(!is.na(str.CvD.res$EntrezGene)),]
row.names(str.CvD.res) <- str.CvD.res$EntrezGene
str.CvD.res$EntrezGene <- NULL
colnames(str.CvD.res)[1] <- "RefSeq_ID"
str.CvD.res <- str.CvD.res[,c((ncol(str.CvD.res)-1),1,(ncol(str.CvD.res)),seq(from=2,to=(ncol(str.CvD.res)-2)))]

#Order by adjusted p value and fold change
str.CvD.res <- str.CvD.res[order(str.CvD.res$padj,-abs(str.CvD.res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(str.CvD.res$padj < 0.1))
length(which(str.CvD.res$pvalue < 0.05))

#Export the results
write.table(transform(as.data.frame(str.CvD.res), entrez_ID = rownames(str.CvD.res))[,c(length(colnames(str.CvD.res))+1,1:length(colnames(str.CvD.res)))],file=paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(str.CvD.res, outliers = T, maxXlim = 2, minXlim = -2, autoScaleAxes = T, maxYlim = 7.5, labels = F, autoScaleLabels = F, maxLabels = 50, log2FC = 0.5)
ggsave(paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
dir.create(paste("str.CvD.differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

#Plot the top 10 differential genes from our results, and save
for (theGene in str.CvD.res[1:10,"Symbol"]){
  ggPlotCountsGR(theGene, intgroup = "GR", subgroup=1:length(colnames(str.CvD.dds)), str.CvD.res, str.CvD.dds, str.CvD.colData)
  ggsave(paste("str.CvD.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#------------------------------------------------------------------------------------------
#RPKM heatmap
nrow(subset(str.CvD.res, abs(log2FoldChange) > 3))
pheatmap(subset(str.CvD.res, abs(log2FoldChange) > 3)[,c(str.CvD.PATIENT_TREATED, str.CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(str.CvD.res, abs(log2FoldChange) > 3)[,"Symbol"],
         annotation_col = as.data.frame(colData(str.CvD.dds)[,c("Group", "Sex", "Response")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.log2FC3.pdf",  sep = "")
)

pheatmap(subset(str.CvD.res, pvalue < 0.05 & abs(log2FoldChange) > 2)[,c(str.CvD.PATIENT_TREATED, str.CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(str.CvD.res, pvalue < 0.05 & abs(log2FoldChange) > 2)[,"Symbol"],
         annotation_col = as.data.frame(colData(str.CvD.dds)[,c("Group", "Sex", "Response")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.pval0.05log2FC2.pdf",  sep = "")
)

pheatmap(subset(str.CvD.res, pvalue<0.05)[,c(str.CvD.PATIENT_TREATED, str.CvD.PATIENT_HIGH)],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = as.data.frame(colData(str.CvD.dds)[,c("Group", "Sex", "Response")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 1,
         cellwidth = 20,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.pvalue0.05.pdf",  sep = "")
)


#Make heatmap of select genes
genes <- paste("^", c("TNF", "IL1B", "IL6", "CXCL8", "IL18", "CCL2", "CD36", "MMP2", "MMP8", "MMP9" , "MARCO" ), "$", sep = "")
str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), str.CvD.res$Symbol)]

pheatmap(subset(str.CvD.res, Symbol %in% str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), str.CvD.res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(str.CvD.res, Symbol %in% str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), str.CvD.res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("str.CvD.differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = "")
)
#And plot them
dir.create(paste("str.CvD.differential_expression/", contrast.name, "/select INF gene counts/", sep = ""), showWarnings = FALSE)
listofgenes <- str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), str.CvD.res$Symbol)]
for (theGene in listofgenes){
  ggPlotCountsGR(theGene, intgroup = "GR", subgroup=1:length(colnames(str.CvD.dds)), str.CvD.res, str.CvD.dds, str.CvD.colData)
  ggsave(paste("str.CvD.differential_expression/", contrast.name, "/select INF gene counts/", theGene, " counts.pdf", sep = ""), width = 10, height = 10)
}

#Also plot these genes for the general groups
dir.create(paste("select INF genes differential_expression/", contrast.name, "/select INF gene counts/", sep = ""), showWarnings = FALSE)
listofgenes <- str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), str.CvD.res$Symbol)]
for (theGene in listofgenes){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("select INF genes differential_expression/", contrast.name, "/select INF gene counts/", theGene, " counts.pdf", sep = ""), width = 10, height = 10)
}

res[which(res$Symbol == "IL6"),]
#Check pathways
get_ontology(res = str.CvD.res, name = paste(contrast.name, "/GSEA", sep = ""), mydir = "str.CvD.differential_expression/", padj.cutoff = 1, log2FC.cutoff = 2)



