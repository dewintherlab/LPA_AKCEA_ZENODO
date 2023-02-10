#Load Packages
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(BiocParallel)
library(genefilter)
library(ggplot2)
library(grid)
library(scatterplot3d)
library(pheatmap)
library(reshape2)
library(org.Hs.eg.db)
library(edgeR)
library(EGSEA)
library(WGCNA)
library(ggrepel)
library(biomaRt)
library(ape)
library(phytools)
allowWGCNAThreads()
register(MulticoreParam(4))
source("functions.R")

#-----------------------------------------------------------------------------------------
#Set up gene ontology objects for pathway analyses
ont_sets <- readList("cp_and_h.gmt.txt")

marts <- listMarts()
mart <- marts[1,1]

human <- useMart(mart, dataset = "hsapiens_gene_ensembl")

#------------------------------------------------------------------------------------------
#Setup heatmap annotation colors
ann_colors <- list(
  Sex = c(F = "darksalmon",M = "cornflowerblue"),
  LPA = c(High = "darkseagreen4", Low = "chocolate1"),
  FlowCell = c(HCLLKDSXX = "firebrick1", HCMKVDSXX = "firebrick4")
)

#------------------------------------------------------------------------------------------
#Load Data
countData <- read.table("Processed_data/LPA.raw_read_quant.table.clean.txt",header=T,row.names=1)
rpkmData <- read.table("Processed_data/LPA.RPKM_quant.table.clean.txt",header=T,row.names=1)
colData <- read.table("Processed_data/metadata.table.txt",header=T,row.names=1, colClasses = "factor") #Added imputed gender based on chrY counts

#Remove outliers S01, S02, S03, S09, S17, S11
rpkmData <- rpkmData[,which(colnames(rpkmData) != "S01")]
countData <- countData[,which(colnames(countData) != "S01")]
colData <- colData[which(row.names(colData) != "S01"),]

rpkmData <- rpkmData[,which(colnames(rpkmData) != "S02")]
countData <- countData[,which(colnames(countData) != "S02")]
colData <- colData[which(row.names(colData) != "S02"),]

rpkmData <- rpkmData[,which(colnames(rpkmData) != "S03")]
countData <- countData[,which(colnames(countData) != "S03")]
colData <- colData[which(row.names(colData) != "S03"),]

rpkmData <- rpkmData[,which(colnames(rpkmData) != "S09")]
countData <- countData[,which(colnames(countData) != "S09")]
colData <- colData[which(row.names(colData) != "S09"),]

rpkmData <- rpkmData[,which(colnames(rpkmData) != "S17")]
countData <- countData[,which(colnames(countData) != "S17")]
colData <- colData[which(row.names(colData) != "S17"),]

rpkmData <- rpkmData[,which(colnames(rpkmData) != "S11")]
countData <- countData[,which(colnames(countData) != "S11")]
colData <- colData[which(row.names(colData) != "S11"),]

#Remove erronous samples P63 and P67 (S02, S06)
#S02 was already removed
rpkmData <- rpkmData[,which(colnames(rpkmData) != "S06")]
countData <- countData[,which(colnames(countData) != "S06")]
colData <- colData[which(row.names(colData) != "S06"),]

#Rectify wrong classification of P73 (S11)
#colData["S11","LPA"] <- "High"

#Set groups
controls <- which(colData$LPA == "Low")
patients <- which(colData$LPA == "High")

#-----------------------------------------------------------------------------------------
#Remove low or high count genes
#Keep median RPKM > 0 in at least one group
keep <- apply(rpkmData[,controls],1,median) > 0 | apply(rpkmData[,patients],1,median) > 0
rpkmData <- rpkmData[keep,]
countData <- countData[keep,]

#Remove median RPKM > 2000
keep <- apply(rpkmData[,c(controls,patients)],1,median) < 2000
rpkmData <- rpkmData[keep,]
countData <- countData[keep,]

#-----------------------------------------------------------------------------------------
#Clean up the tables
#Make all counts integers
countData <- as.data.frame(apply(countData, 1:2, as.integer))

#Remove NAs
indx <- apply(countData, 1, function(x) any(is.na(x)))
countData <- countData[-!indx,]

indx <- apply(rpkmData, 1, function(x) any(is.na(x)))
rpkmData <- rpkmData[-!indx,]

#-----------------------------------------------------------------------------------------
#Annotate RPKM table
egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(row.names(rpkmData),egREFSEQ$accession)
rpkmData$EntrezGene <- egREFSEQ$gene_id[m]

egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(rpkmData$EntrezGene, egSYMBOL$gene_id)
rpkmData$Symbol <- egSYMBOL$symbol[m]

egCHR <- toTable(org.Hs.egCHR)
m <- match(rpkmData$EntrezGene, egCHR$gene_id)
rpkmData$Chr <- egCHR$chromosome[m]

#-----------------------------------------------------------------------------------------
#Remove all but the highest expressed transcript per gene
o <- order(rowSums(rpkmData[,1:(length(colnames(rpkmData))-3)]), decreasing=TRUE)
countData <- countData[o,]
rpkmData <- rpkmData[o,]
d <- duplicated(rpkmData$Symbol)
countData <- countData[!d,]
rpkmData <- rpkmData[!d,]

#Remove ribosomal, mitochondrial, and MIR genes
keep <- grep("^MT-|^RPL|^RPS|^MIR|^SNOR|^SCAR", rpkmData$Symbol, invert = T)
rpkmData <- rpkmData[keep,]
countData <- countData[keep,]

#Remove chrom X and Y genes
#keep <- grep("X|Y", rpkmData$Chr, invert = T)
#rpkmData <- rpkmData[keep,]
#countData <- countData[keep,]

#------------------------------------------------------------------------------------------
#Setup data
#dds <- DESeqDataSetFromMatrix(countData,colData,design= ~Sex + LPA)
dds <- DESeqDataSetFromMatrix(countData,colData,design= ~LPA)
dds$LPA <- relevel(dds$LPA, "Low")
dds <- DESeq(dds,parallel = T)

#Create a directory for the results
dir.create("differential_expression", showWarnings = FALSE)
dir.create("differential_expression_keepsex", showWarnings = FALSE)
dir.create("differential_expression_keepsex_no_outliers_P73_corrected", showWarnings = FALSE)
dir.create("differential_expression_keepsex_no_outliers_P73_gone", showWarnings = FALSE)

#-----------------------------------------------------------------------------------------
#Call differntial expression
res <- results(dds,parallel= T, contrast = c("LPA", "High", "Low"))

#Visualise
pdf("differential_expression_keepsex_no_outliers_P73_gone/LPA_MA_Plot.pdf")
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

#Export the results
write.table(as.data.frame(res),file="differential_expression_keepsex_no_outliers_P73_gone/LPA_diff_genes.DESeq2.txt",quote=F,sep="\t")

#Define new column numbers for groups
controls <- which(colData$LPA == "Low") + 9
patients <- which(colData$LPA == "High") + 9


#------------------------------------------------------------------------------------------
#Make plots
#Setup data for plotting
rld <-rlog(dds)
topVarGenes<-head(order(-rowVars(assay(rld))),1000)
mat <- assay(rld)[topVarGenes,]

#Distance plot
distsRL <- dist(t(assay(rld)))
d.mat <- as.matrix(distsRL)
rownames(d.mat) <- colnames(d.mat) <- paste(colnames(mat),with(colData(dds), paste(LPA, Sex, sep=" : ")),sep = " : ")
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("differential_expression_keepsex_no_outliers_P73_gone/LPA_distance_plot.pdf")
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
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = F,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/LPA_kmeans_heatmap.pdf"
)

#2D PCA plots
PCA.plot(rld, intgroup = "Sex")
ggsave("differential_expression_keepsex_no_outliers_P73_gone/LPA_PCA_plot_Sex.pdf")

PCA.plot(rld, intgroup = "LPA")
ggsave("differential_expression_keepsex_no_outliers_P73_gone/LPA_PCA_plot_LPA.pdf")

#Volcano plot
volcano_plot(res, outliers = T, maxXlim = 2.5, minXlim = -2.5, autoScaleAxes = F, maxYlim = 3, labels = F)
ggsave("differential_expression_keepsex_no_outliers_P73_gone/LPA_volcano.pdf")

#RPKM heatmap
pheatmap(subset(res, padj<0.1 & abs(log2FoldChange)>0.5)[,c(controls,patients)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, padj<0.1 & abs(log2FoldChange)>0.5)[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 2,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/LPA_RPKM_clustered_heatmap_p0.1_FC0.5.pdf"
)

#------------------------------------------------------------------------------------------
#Check counts for some genes
#Plot the top 10 differential genes from our results, and save
dir.create("differential_expression_keepsex_no_outliers_P73_gone/gene counts", showWarnings = FALSE)

for (theGene in res[1:10,"Symbol"]){
  ggPlotCounts(theGene, intgroup = "LPA")
  ggsave(paste("differential_expression_keepsex_no_outliers_P73_gone/gene counts/", theGene, " counts.pdf", sep = ""))
}

ggPlotCounts("NTN1", intgroup = "LPA")
ggsave("differential_expression/gene counts_keepsex/NTN1 counts.pdf")
ggPlotCounts("NTN1", intgroup = "FlowCell")

ggPlotCounts("UNC5B", intgroup = "LPA")
ggsave("differential_expression/gene counts_keepsex/UNC5B counts.pdf")
ggPlotCounts("UNC5B", intgroup = "FlowCell")


pheatmap(res[grep("UNC5B$|NTN1", res$Symbol),c(patients,controls)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = res[grep("UNC5B$|NTN1", res$Symbol),"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/NTN1_UNC5B_RPKM_clustered_heatmap.pdf"
)

#Retrieve and plot GO terms and pathways
get_ontology(res, name = "LPA")

#Full EGSEA report
#Extract genes from res with under cut off (default 0.1)
genes <- data.frame(entrezGene = row.names(subset(res, padj < 0.1)), log2FoldChange = subset(res, padj < 0.1)$log2FoldChange)

#Make symbols object
gene_symbols <- data.frame(entrezGene = genes$entrezGene, Symbol = subset(res, padj < 0.1)$Symbol)

#Build annotation index
gs.annots <- buildIdx(entrezIDs = genes$entrezGene,
                        species = "human",
                   msigdb.gsets = "all",
                     gsdb.gsets = "all",
                        go.part = T,
                   kegg.updated = T
)

#Run EGSEA
egsea.ora(geneIDs = genes$entrezGene,
            logFC = genes$log2FoldChange,
            title = "LPA", 
         universe = rpkmData$EntrezGene,
        gs.annots = gs.annots,
       symbolsMap = gene_symbols,
      display.top = 10,
          sort.by = "p.adj",
       report.dir = "differential_expression_keepsex_no_outliers_P73_gone/EGSEA",
         kegg.dir = "differential_expression_keepsex_no_outliers_P73_gone/EGSEA/kegg-dir",
      num.threads = 4,
      interactive = T,
           report = T,
          verbose = F
)


#Make heatmap of select genes
genes <- paste("^", c("CCR7", "SELL", "ITGAM", "ITGAX", "ITGB1", "CD36", "SRA1", "CD163","CD200R1", "MRC1", "IL1B", "IL6", "TNF", "IL10"), "$", sep = "")
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,c(controls,patients)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/LPA_RPKM_clustered_heatmap_select_genes.pdf")


#Make heatmap of plasma proteo genes
genes <- paste("^", c("OSM", "CD6", "MARCO","TNFSH10A", "CASP8","EIF4EBP1","THBD","CXCL11",
                      "HB-EGF", "CCL8", "CD5","CD84","CEACAM8", "CD40LG"),
               "$", sep = "")
genes <- unique(genes)
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,c(controls,patients)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/LPA_RPKM_clustered_heatmap_select_genes_from_plasma_proteo.pdf")

#Make heatmap of cell lysate proteo genes
genes <- paste("^", c("MARCO","HB-EGF", "SPON2",  "CX3CL1", "CCL3" , "IL10-RA","CCl4",
                      "TGFb-1", "CCl20" , "CXCL11","CCL23" ,"CCL2" , "CST5"),
               "$", sep = "")
genes <- unique(genes)
res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)]

pheatmap(subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,c(controls,patients)],
         scale = "row",
         show_rownames = T,
         show_colnames = T,
         labels_row = subset(res, Symbol %in% res$Symbol[grep(paste(genes, collapse = "|"), res$Symbol)])[,"Symbol"],
         annotation_col = as.data.frame(colData(dds)[,c("LPA", "Sex")]),
         annotation_colors = ann_colors,
         drop_levels = T,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row = 7,
         cutree_cols = 1,
         cutree_rows = 1,
         cluster_rows = F,
         filename = "differential_expression_keepsex_no_outliers_P73_gone/LPA_RPKM_clustered_heatmap_select_genes_from_cell_lysate_proteo.pdf")
