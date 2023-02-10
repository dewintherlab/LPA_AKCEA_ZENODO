#Setup dir
contrast.name <- "ALL_DIFFS"
dir.create(paste("differential_expression/", contrast.name, sep = ""), showWarnings = FALSE)

#Flatten lists and retrieve RPKM object
all.diffs.uniq <- unique(unlist(all.diffs))
all.diffs.RPKM <- rpkmData[rpkmData$Symbol %in% all.diffs.uniq,]
RPKM.cols <- 1:(length(colnames(all.diffs.RPKM)) - 3 )

#Define column numbers for RPKM object
HEALTHY_BASELINE <- which(colData$Group == "Healthy_baseline_LPA")
HEALTHY_HIGH     <- which(colData$Group == "Healthy_high_LPA")
PATIENT_HIGH     <- which(colData$Group == "Patient_high_LPA")
PATIENT_TREATED  <- which(colData$Group == "Patient_treated_LPA")

#RPKM Heatmap
pheatmap(all.diffs.RPKM[,RPKM.cols],
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = data.frame(Flowcell=colData$Flowcell, Sex=colData$Sex, Group=colData$Group, row.names = row.names(colData)),
         annotation_colors = ann_colors,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = F,
         kmeans_k = 6,
         cellwidth = 10,
         cellheight = 40,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".kmeans_heatmap.pdf",  sep = "")
)

#RPKM Heatmap of most variable genes
topVarGenes <- head(order(-rowVars(all.diffs.RPKM[,RPKM.cols])),1000)
pheatmap(all.diffs.RPKM[topVarGenes,RPKM.cols],
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = data.frame(Flowcell=colData$Flowcell, Sex=colData$Sex, Group=colData$Group, row.names = row.names(colData)),
         annotation_colors = ann_colors,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = F,
         cellwidth = 10,
         cellheight = 1,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, "top_1000_var.RPKM_heatmap.pdf",  sep = "")
)

#Euclidean Distance plot
distsRL <- dist(t(all.diffs.RPKM[,RPKM.cols]))
d.mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pheatmap(d.mat,
         cluster_rows = hc, 
         cluster_cols = hc,
         color = rev(hmcol),
         cellwidth = 10,
         cellheight = 10,
         border_color = NA,
         annotation_row = data.frame(Flowcell=colData$Flowcell, Sex=colData$Sex, Group=colData$Group, row.names = row.names(colData)),
         annotation_col = data.frame(Flowcell=colData$Flowcell, Sex=colData$Sex, Group=colData$Group, row.names = row.names(colData)),
         annotation_colors = ann_colors,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".distance_plot.pdf",  sep = "")
        )


#-------------------------------------------------------------------------------------------------------------------------
#Read top genes from AKCEA
AKCEA.treatment_affected_genes <- scan(file = "BASELINE_v_TREATMENT.significant_genes_padj0.1_log2FC0.5.txt", what = "list")
topVarGenes <- rpkmData$Symbol %in% AKCEA.treatment_affected_genes

pheatmap(rpkmData[topVarGenes,RPKM.cols],
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = data.frame(Flowcell=colData$Flowcell, Sex=colData$Sex, Group=colData$Group, row.names = row.names(colData)),
         annotation_colors = ann_colors,
         annotation_names_row = F,
         show_rownames = T,
         show_colnames = F,
         cellwidth = 10,
         cellheight = 10,
         filename = paste("differential_expression/", contrast.name, "/", contrast.name, ".AKCEA.genes.RPKM_heatmap.pdf",  sep = "")
)

#Average RPKM heatmap
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED], 1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Patient Treated","Healthy Baseline","Healthy High LPA")
topdiff.res <- as.matrix(topdiff.res)

topdiff.res <- topdiff.res[which(apply(topdiff.res,1,max) < 100),] #Only print RPKMs < 100 for improved visibility


pheatmap(topdiff.res[sort(topdiff.res[,1], index.return=T)$ix,],
         cluster_rows = F,
         cluster_cols = F,
         labels_col =  colnames(topdiff.res),
         labels_row = rpkmData[topVarGenes,][sort(topdiff.res[,1], index.return=T)$ix,"Symbol"],
         cellwidth = 25,
         cellheight = 12,
         border_color = NA,
         angle_col = 45, display_numbers = T,
         filename =  paste("differential_expression/", contrast.name, "/", contrast.name, ".AKCEA.genes.Median_RPKM_heatmap.pdf", sep = ""),
         col=colorRampPalette(c("white","chartreuse2","salmon","firebrick2"))(nrow(topdiff.res))
)

m <- melt(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Patient Treated"]),])
colnames(m) <- c("refSeq_ID","Group","RPKM")
m$Symbol <- rpkmData[topVarGenes,][rownames(rpkmData[topVarGenes,]) %in% m$refSeq_ID,"Symbol"]

ggplot(m,aes(x = Group, y = RPKM, col=Symbol, group = Symbol)) +
  geom_point(size=3) +
  geom_line() +
  scale_y_log10() +
  ylab("log10(RPKM)") +
  ggtitle("Median RPKM of differentially expressed IFN gamma signaling genes") +
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
ggsave(paste("differential_expression/", contrast.name, "/", contrast.name, ".AKCEA.genes.Median_RPKM_lowered_plot.pdf", sep = ""))



#Average plots
#Get top 1000 most varying genes
topVarGenes <- head(order(-rowVars(all.diffs.RPKM[,RPKM.cols])),1000)

#Take median per gene
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED], 1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Patient Treated","Healthy Baseline","Healthy High LPA")
topdiff.res <- as.matrix(topdiff.res)

#Genes down after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Patient Treated"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down after treatment") +
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

#Genes up after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] < topdiff.res[,"Patient Treated"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up after treatment") +
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

#Genes low in baseline
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] < topdiff.res[,"Healthy High LPA"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes low at baseline") +
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

#Genes high in baseline
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] > topdiff.res[,"Healthy High LPA"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes high at baseline") +
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


#Genes down 1.5x after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] / topdiff.res[,"Patient Treated"] >= 1.5),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down 1.5x after treatment") +
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


#Genes up 1.4x after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient Treated"] / topdiff.res[,"Patient High LPA"] >= 1.4),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up 1.4x after treatment") +
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

#Genes up 2x in patients
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] / topdiff.res[,"Healthy Baseline"] >= 2),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up 2x in patients") +
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

#Genes down 2x in patients
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] / topdiff.res[,"Patient High LPA"] >= 2),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down 2x in patients") +
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



#Read top genes from AKCEA
topVarGenes <- rpkmData$Symbol %in% AKCEA.treatment_affected_genes

#Take median per gene
topdiff.res     <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_HIGH],    1,median))
topdiff.res[,2] <- as.data.frame(apply(rpkmData[topVarGenes,PATIENT_TREATED], 1,median))
topdiff.res[,3] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_BASELINE],1,median))
topdiff.res[,4] <- as.data.frame(apply(rpkmData[topVarGenes,HEALTHY_HIGH],    1,median))
colnames(topdiff.res) <- c("Patient High LPA","Patient Treated","Healthy Baseline","Healthy High LPA")
topdiff.res <- as.matrix(topdiff.res)

#Genes down after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] > topdiff.res[,"Patient Treated"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down after treatment") +
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

#Genes up after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] < topdiff.res[,"Patient Treated"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up after treatment") +
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

#Genes low in baseline
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] < topdiff.res[,"Healthy High LPA"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes low at baseline") +
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

#Genes high in baseline
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] > topdiff.res[,"Healthy High LPA"]),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes high at baseline") +
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


#Genes down 1.5x after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] / topdiff.res[,"Patient Treated"] >= 1.5),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down 1.5x after treatment") +
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


#Genes up 1.4x after treatment
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient Treated"] / topdiff.res[,"Patient High LPA"] >= 1.4),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up 1.4x after treatment") +
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

#Genes up 2x in patients
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Patient High LPA"] / topdiff.res[,"Healthy Baseline"] >= 2),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes up 2x in patients") +
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

#Genes down 2x in patients
#Take stats per group
topdiff.res.summ <- apply(topdiff.res[which(topdiff.res[,"Healthy Baseline"] / topdiff.res[,"Patient High LPA"] >= 2),],2,summary)

#Prep data for plot
m     <- melt(topdiff.res.summ[c(3),])     #Get median
m$q25 <- melt(topdiff.res.summ[c(2),])[,1] #Get 25th percentile
m$q75 <- melt(topdiff.res.summ[c(5),])[,1] #Get 75th percentile
m$group <- factor(row.names(m), levels=row.names(m)) #Set groups
colnames(m) <- c("RPKM","q25","q75","group")

ggplot(m, aes(x = group, y = RPKM, ymin = q25, ymax = q75)) +
  geom_ribbon(group = "RPKM", fill = "dodgerblue", alpha = 0.25) +
  geom_line(group="RPKM", color = "dodgerblue", size = 1.5) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,2.1)) +
  ggtitle("Average RPKM of genes down 2x in patients") +
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


#------------------------------------------------------------------------------------
#Make barplot of DE genes in the various groups
DEs <- data.frame(AvB=c(43, 52),
                  AvC=c(769, 517),
                  CvD=c(359, 216),
                  row.names = c("Up", "Down"))

m <- melt(DEs)
m$Direction <- c("Up", "Down")
colnames(m) <- c("Group", "DE", "Direction")
m$Direction <- factor(m$Direction, levels = c("Up", "Down"))

ggplot(m, aes(x = Group, y = DE, fill = Direction)) +
  geom_bar(stat = "identity",  position = "dodge") +
  scale_fill_manual("Direction", values = c(Up = "coral3", Down = "cornflowerblue")) +
  scale_x_discrete(labels = c("Healthy controls - High vs. Low LP(a)", 
                             "CVD Patients high LP(a) vs. Healthy Controls - low LP(a)", 
                             "CVD Patients high LP(a) vs. CVD Patients ASO treated low LP(a)")) +
  theme_light() +
  ylab("# of DE genes") +
  scale_y_continuous(limits = c(0, 800), breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800)) +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        aspect.ratio = 2/1)
ggsave("number of DE genes per comparison.pdf")
