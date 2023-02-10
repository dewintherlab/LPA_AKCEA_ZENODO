#Miscellaneous graphs
#------------------------------------------------------------------
dir.create("misc_graphs", showWarnings = FALSE)
RPKM.cols <- seq(10,length(colnames(res)))

#Setup data
r <- res[,RPKM.cols]
row.names(r) <- res$Symbol
ra <- data.frame("Healthy_Baseline" = apply(r[,rownames(colData[which(colData$Group == "Healthy_baseline_LPA"),])],1,mean), 
                 "Healthy_High_LPA" = apply(r[,rownames(colData[which(colData$Group == "Healthy_high_LPA"),])],1,mean),
                 "Patient_High_LPA" = apply(r[,rownames(colData[which(colData$Group == "Patient_high_LPA"),])],1,mean),
              "Patient_Treated_LPA" = apply(r[,rownames(colData[which(colData$Group == "Patient_treated_LPA"),])],1,mean)
              )

#PPathway heatmaps
#Migration pathway
MIG <- scan("KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION.txt", what = "list")
rs  <-  r[rownames(r)  %in% MIG,]
ras <- ra[rownames(ra) %in% MIG,]
ras <- ras[apply(ras,1,sum) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION.RPKM.heatmap.pdf"
         )

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION.average_RPKM.heatmap.pdf"
         )

#Monocyte pathway
MON <- scan("BIOCARTA_MONOCYTE_PATHWAY.txt", what = "list")
rs  <-  r[rownames(r)  %in% MON,]
ras <- ra[rownames(ra) %in% MON,]
ras <- ras[apply(ras,1,sum) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/BIOCARTA_MONOCYTE_PATHWAY.RPKM.heatmap.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/BIOCARTA_MONOCYTE_PATHWAY.average_RPKM.heatmap.pdf"
)

#Adhesion pathway
ADH <- scan("KEGG_CELL_ADHESION_MOLECULES_CAMS.txt", what = "list")
rs  <-  r[rownames(r)  %in% ADH,]
ras <- ra[rownames(ra) %in% ADH,]
ras <- ras[apply(ras,1,sum) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/KEGG_CELL_ADHESION_MOLECULES_CAMS.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_CELL_ADHESION_MOLECULES_CAMS.average_RPKM.heatmap.pdf"
)

#Cyotkine pathway
CYT <- scan("BIOCARTA_CYTOKINE_PATHWAY.txt", what = "list")
rs  <-  r[rownames(r)  %in% CYT,]
ras <- ra[rownames(ra) %in% CYT,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/BIOCARTA_CYTOKINE_PATHWAY.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/BIOCARTA_CYTOKINE_PATHWAY.average_RPKM.heatmap.pdf"
)

#Inflammation pathway
INFL <- scan("BIOCARTA_INFLAM_PATHWAY.txt", what = "list")
rs  <-  r[rownames(r)  %in% INFL,]
ras <- ra[rownames(ra) %in% INFL,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/BIOCARTA_INFLAM_PATHWAY.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/BIOCARTA_INFLAM_PATHWAY.average_RPKM.heatmap.pdf"
)

#TLR pathway
TLR <- scan("KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY.txt", what = "list")
rs  <-  r[rownames(r)  %in% TLR,]
ras <- ra[rownames(ra) %in% TLR,]
ras <- ras[apply(ras,1,max) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY.average_RPKM.heatmap.pdf"
)


#IFN singaling
IFN <- scan("REACTOME_INTERFERON_SIGNALING.txt", what = "list")
rs  <-  r[rownames(r)  %in% IFN,]
ras <- ra[rownames(ra) %in% IFN,]
ras <- ras[apply(ras,1,max) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/REACTOME_INTERFERON_SIGNALING.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/REACTOME_INTERFERON_SIGNALING.average_RPKM.heatmap.pdf"
)


#Adjust cut-off
ras <- ras[apply(ras,1,max) > 5,]
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/REACTOME_INTERFERON_SIGNALING.OVER5_average_RPKM.heatmap.pdf"
)


#KEGG chemokine singaling
KCSP <- scan("KEGG_CHEMOKINE_SIGNALING_PATHWAY.txt", what = "list")
rs  <-  r[rownames(r)  %in% KCSP,]
ras <- ra[rownames(ra) %in% KCSP,]
ras <- ras[apply(ras,1,max) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/KEGG_CHEMOKINE_SIGNALING_PATHWAY.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_CHEMOKINE_SIGNALING_PATHWAY.average_RPKM.heatmap.pdf"
)


#Adjust cut-off
ras <- ras[apply(ras,1,max) > 5,]
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_CHEMOKINE_SIGNALING_PATHWAY.OVER5_average_RPKM.heatmap.pdf"
)

#KEGG cytokine interaction singaling
KCCRI <- scan("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.txt", what = "list")
rs  <-  r[rownames(r)  %in% KCCRI,]
ras <- ra[rownames(ra) %in% KCCRI,]
ras <- ras[apply(ras,1,max) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.average_RPKM.heatmap.pdf"
)


#Adjust cut-off
ras <- ras[apply(ras,1,max) > 5,]
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.OVER5_average_RPKM.heatmap.pdf"
)



#REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES singaling
RCRBC <- scan("REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES.txt", what = "list")
rs  <-  r[rownames(r)  %in% RCRBC,]
ras <- ra[rownames(ra) %in% RCRBC,]
ras <- ras[apply(ras,1,max) > 1,]
ras <- ras[ras[,1] < 100,]

#Make RPKM heatmap
pheatmap(rs, scale = "row",
         show_rownames = T,
         show_colnames = F, 
         cellheight = 10, 
         annotation_col = colData[,"Group",drop =F],
         filename = "misc_graphs/REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES.pdf"
)

#Make gradient heatmap of averages
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES.average_RPKM.heatmap.pdf"
)


#Adjust cut-off
ras <- ras[apply(ras,1,max) > 5,]
pheatmap(ras[sort(ras[,1], index.return=T)$ix,], 
         cluster_rows = F, 
         cluster_cols = F,
         color = colorRampPalette(c("green","yellow","goldenrod","darkorange","red","darkred"))(nrow(ras)),
         cellwidth = 30,
         cellheight = 10,
         border_color = NA,
         display_numbers = T,
         filename = "misc_graphs/REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES.OVER5_average_RPKM.heatmap.pdf"
)

dir.create("CYT_genes", showWarnings = FALSE)
rs <- rownames(r)[rownames(r)  %in% CYT]

for (gene in rs){
  ggPlotCounts(theGene = gene,intgroup = "Group")
  ggsave(file = paste("CYT_genes/",gene, "_norm_count.pdf", sep = ""))
}

ras <- ra[rownames(ra) %in% CYT,]
ras <- apply(ras,1, mean)
CYT.df <- data.frame(Expressed=names(ras[ras >1]),
                     Lowly_Expressed=names(ras[ras > 0.1 & ras < 1]),
                     Barely_Detected=names(ras[ras < 0.1]),
                     Not_Expressed=CYT[!(CYT %in% rs)]
                     )
write.table(CYT.df, file = paste("CYT_genes/", gene, "CYT_GENES_Expression.txt", sep = ""), quote = F, row.names = F, sep = "\t")
