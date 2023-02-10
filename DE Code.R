#Call DE genes for all relevant contrasts
dir.create("differential_expression", showWarnings = FALSE)

#Define new column numbers for results object
ALL.HEALTHY_BASELINE <- which(colData$Group == "Healthy_baseline_LPA") + 9
ALL.HEALTHY_HIGH     <- which(colData$Group == "Healthy_high_LPA")     + 9
ALL.PATIENT_HIGH     <- which(colData$Group == "Patient_high_LPA")     + 9
ALL.PATIENT_TREATED  <- which(colData$Group == "Patient_treated_LPA")  + 9

RPKM.cols <- c(10:(length(colnames(countData)) + 9))

#------------------------------------------------------------------------------------------
#Call differential expression for healthy baseline v. high
contrast.name <- "HEALTHY_BASELINE_v_HIGH"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Healthy_high_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = T, maxXlim = 25, minXlim = -25, autoScaleAxes = F, maxYlim = 3, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

ggPlotCounts("NTN1", intgroup = "Group")
ggsave(paste("differential_expression/", contrast.name, "/gene counts/", "NTN1", " counts.pdf", sep = ""))

ggPlotCounts("UNC5B", intgroup = "Group")
ggsave(paste("differential_expression/", contrast.name, "/gene counts/", "UNC5B", " counts.pdf", sep = ""))

pheatmap(res[grep("^UNC5B$|^NTN1$", res$Symbol),c(RPKM.cols)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = res[grep("UNC5B$|NTN1", res$Symbol),"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         filename = paste("differential_expression/", contrast.name, "/UNC5B_NTN1.heatmap.pdf", sep = "")
)


#Make heatmap of select genes
genes <- paste("^", c("CCR7", "SELL", "ITGAM", "ITGAX", "ITGB1", "CD36", "SRA1", "CD163","CD200R1", "MRC1", "IL1B", "IL6", "TNF", "IL10"), "$", sep = "")
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("differential_expression/", contrast.name, "/select_genes.heatmap.pdf", sep = ""))


#Make heatmap of plasma proteo genes
genes <- paste("^", c("OSM", "CD6", "MARCO","TNFSH10A", "CASP8","EIF4EBP1","THBD","CXCL11",
                      "HB-EGF", "CCL8", "CD5","CD84","CEACAM8", "CD40LG"),
               "$", sep = "")
genes <- unique(genes)
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("differential_expression/", contrast.name, "/plasma_proteo_genes.heatmap.pdf", sep = ""))

#Make heatmap of cell lysate proteo genes
genes <- paste("^", c("MARCO","HB-EGF", "SPON2",  "CX3CL1", "CCL3" , "IL10-RA","CCl4",
                      "TGFb-1", "CCl20" , "CXCL11","CCL23" ,"CCL2" , "CST5"),
               "$", sep = "")
genes <- unique(genes)
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,RPKM.cols],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = paste("differential_expression/", contrast.name, "/cell_lysate_proteo_genes.heatmap.pdf", sep = ""))


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient high v. treated
contrast.name <- "PATIENT_TREATED_v_HIGH"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_treated_LPA", "Patient_high_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 2.5, minXlim = -2.5, autoScaleAxes = F, maxYlim = 3, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

ggPlotCounts("ATF7", intgroup = "Group")

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.05 & abs(log2FoldChange)>0.1)[,c(ALL.PATIENT_HIGH,ALL.PATIENT_TREATED)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.05 & abs(log2FoldChange)>0.1)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.05.log2FC0.1.pdf",  sep = "")
)

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient high v. healthy baseline
contrast.name <- "PATIENT_HIGH_v_HEALTHY_BASELINE"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_high_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 5, minXlim = -5, autoScaleAxes = F, maxYlim = 4, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.01 & abs(log2FoldChange)>1)[,c(ALL.PATIENT_HIGH,ALL.HEALTHY_BASELINE)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.01 & abs(log2FoldChange)>1)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.01.log2FC1.pdf",  sep = "")
)

#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient treated v. healthy baseline
contrast.name <- "PATIENT_TREATED_v_HEALTHY_BASELINE"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_treated_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 2.5, minXlim = -2.5, autoScaleAxes = F, maxYlim = 3, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.1)[,c(ALL.PATIENT_TREATED, ALL.HEALTHY_BASELINE)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.1)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.pdf",  sep = "")
)


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient treated v. healthy high
contrast.name <- "PATIENT_TREATED_v_HEALTHY_HIGH"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_treated_LPA", "Healthy_high_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 5, minXlim = -5, autoScaleAxes = F, maxYlim = 4, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.1 & abs(log2FoldChange)>1)[,c(ALL.PATIENT_TREATED, ALL.HEALTHY_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.1 & abs(log2FoldChange)>1)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.log2FC1.pdf",  sep = "")
)


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient treated v. healthy baseline
contrast.name <- "PATIENT_TREATED_v_HEALTHY_BASELINE"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_treated_LPA", "Healthy_baseline_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 2.5, minXlim = -2.5, autoScaleAxes = F, maxYlim = 3, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.1)[,c(ALL.PATIENT_TREATED, ALL.HEALTHY_BASELINE)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.1)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.1.pdf",  sep = "")
)


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Call differential expression for patient high v. healthy high
contrast.name <- "PATIENT_HIGH_v_HEALTHY_HIGH"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)
imagemat(model.matrix(design, colData))
colnames(model.matrix(design, colData))
resultsNames(dds)

res <- results(dds, parallel= T, contrast = c("Group", "Patient_high_LPA", "Healthy_high_LPA"))

#Visualise
pdf(paste("differential_expression/", contrast.name, "/", contrast.name, ".MA_Plot.pdf", sep = ""))
DESeq2::plotMA(res, ylim = c(-5, 5))
dev.off()

#Add RPKM values to the results table
res <- merge(as.data.frame(res),rpkmData,by=0)

#Clean up results table
res <- res[which(!is.na(res$EntrezGene)),]
row.names(res) <- res$EntrezGene
res$EntrezGene <- NULL
colnames(res)[1] <- "RefSeq_ID"
res <- res[,c((ncol(res)-1),1,(ncol(res)),seq(from=2,to=(ncol(res)-2)))]

#Order by adjusted p value and fold change
res <- res[order(res$padj,-abs(res$log2FoldChange)),]

#How many diff. genes have we got?
length(which(res$padj < 0.1))
length(which(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5))

#Save the list for later
#all.diffs[contrast.name] <- list(subset(res, padj < 0.1)$Symbol)

#Export the results
write.table(transform(as.data.frame(res), entrez_ID = rownames(res))[,c(length(colnames(res))+1,1:length(colnames(res)))],file=paste("differential_expression/", contrast.name, "/", contrast.name, ".diff_genes.DESeq2.txt",  sep = ""), row.names = F, quote=F, sep="\t")

#Volcano plot
volcano_plot(res, outliers = F, maxXlim = 5, minXlim = -5, autoScaleAxes = F, maxYlim = 7.5, labels = F)
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".volcano_plot.pdf",  sep = ""))


#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create(paste("differential_expression/", contrast.name, "/gene counts/", sep = ""), showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "Group")
  ggsave(paste("differential_expression/", contrast.name, "/gene counts/", theGene, " counts.pdf", sep = ""))
}

get_ontology(res = res, name = paste(contrast.name, "/GSEA", sep = ""))

#RPKM heatmap
pheatmap(subset(res, padj<0.05 & abs(log2FoldChange)>1.5)[,c(ALL.PATIENT_HIGH, ALL.HEALTHY_HIGH)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.05 & abs(log2FoldChange)>1.5)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("Group", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".RPKM_heatmap.padj0.05.log2FC1.5.pdf",  sep = "")
)
