#Call DE genes for all relevant contrasts
dir.create("ultra.str.CvD.differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
ultra.str.CvD.RESPONSE_HIGH <- which(ultra.str.CvD.colData$Response == "High") + 9
ultra.str.CvD.RESPONSE_LOW     <- which(ultra.str.CvD.colData$Response == "PLow")     + 9

RPKM.cols <- c(10:(length(colnames(ultra.str.CvD.countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "PATIENT_HIGH_RESPONDER_v_PATIENT_LOW_RESPONDER"
design  <- ~Response
dir.create(paste("ultra.str.CvD.differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, ultra.str.CvD.colData))
colnames(model.matrix(design, ultra.str.CvD.colData))
resultsNames(ultra.str.CvD.dds)

ultra.str.CvD.res <- results(ultra.str.CvD.dds, parallel= T, contrast = c("Response", "High", "Low"))

#Visualise
pdf(paste("ultra.str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(ultra.str.CvD.res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
ultra.str.CvD.res <- merge(as.data.frame(ultra.str.CvD.res),ultra.str.CvD.rpkmData,by=0)

#Clean up results table
ultra.str.CvD.res <- ultra.str.CvD.res[which(!is.na(ultra.str.CvD.res$EntrezGene)),]
row.names(ultra.str.CvD.res) <- ultra.str.CvD.res$EntrezGene
ultra.str.CvD.res$EntrezGene <- NULL
colnames(ultra.str.CvD.res)[1] <- "RefSeq_ID"
ultra.str.CvD.res <- ultra.str.CvD.res[,c((ncol(ultra.str.CvD.res)-1),1,(ncol(ultra.str.CvD.res)),seq(from=2,to=(ncol(ultra.str.CvD.res)-2)))]

#Order by adjusted p value and fold change
ultra.str.CvD.res <- ultra.str.CvD.res[order(ultra.str.CvD.res$padj,-abs(ultra.str.CvD.res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(ultra.str.CvD.res$padj < 0.1))
length(which(ultra.str.CvD.res$pvalue < 0.05))

#Export the results
write.table(transform(as.data.frame(ultra.str.CvD.res), entrez_ID = rownames(ultra.str.CvD.res))[,c(length(colnames(ultra.str.CvD.res))+1,1:length(colnames(ultra.str.CvD.res)))],file=paste("ultra.str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(ultra.str.CvD.res, outliers = T, maxXlim = 6, minXlim = -6, autoScaleAxes = F, maxYlim = 2.5, labels = T, autoScaleLabels = F, maxLabels = 50, log2FC = 0.5)
ggsave(paste("ultra.str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".labeled.volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
dir.create(paste("ultra.str.CvD.differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

#Plot the top 10 differential genes from our results, and save
for (theGene in ultra.str.CvD.res[1:10,"Symbol"]){
  ggPlotCountsRes(theGene, intgroup = "GR", subgroup=1:length(colnames(ultra.str.CvD.dds)), ultra.str.CvD.res, ultra.str.CvD.dds, ultra.str.CvD.colData)
  ggsave(paste("ultra.str.CvD.differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

#------------------------------------------------------------------------------------------
#RPKM heatmap
nrow(subset(ultra.str.CvD.res, abs(log2FoldChange) > 3))
pheatmap(subset(ultra.str.CvD.res, abs(log2FoldChange) > 3)[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(ultra.str.CvD.res, abs(log2FoldChange) > 3)[,"Symbol"],
         annotation_col = as.data.frame(colData(ultra.str.CvD.dds)[,c("Sex", "Response")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("ultra.str.CvD.differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.log2FC3.pdf",  sep = "")
)

#Make heatmap of select genes
genes <- paste("^", c("TNF", "IL1B", "IL6", "CXCL8", "IL18", "CCL2", "CD36", "MMP2", "MMP8", "MMP9" , "MARCO" ), "$", sep = "")
ultra.str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), ultra.str.CvD.res$Symbol)]

pheatmap(subset(ultra.str.CvD.res, Symbol %in% ultra.str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), ultra.str.CvD.res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(ultra.str.CvD.res, Symbol %in% ultra.str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), ultra.str.CvD.res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("ultra.str.CvD.differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = "")
)
#And plot them
dir.create(paste("ultra.str.CvD.differential_expression/", contrast.name, "/select INF gene counts/", sep = ""), showWarnings = FALSE)
listofgenes <- ultra.str.CvD.res$Symbol[grep(paste(genes, collapse = "|"), ultra.str.CvD.res$Symbol)]
for (theGene in listofgenes){
  ggPlotCountsRes(theGene, intgroup = "GR", subgroup=1:length(colnames(ultra.str.CvD.dds)), ultra.str.CvD.res, ultra.str.CvD.dds, ultra.str.CvD.colData)
  ggsave(paste("ultra.str.CvD.differential_expression/", contrast.name, "/select INF gene counts/", theGene, " counts.pdf", sep = ""))
}


#Check pathways
get_ontology(res = ultra.str.CvD.res, name = paste(contrast.name, "/GSEA", sep = ""), mydir = "ultra.str.CvD.differential_expression/", padj.cutoff = 1, log2FC.cutoff = 2)



