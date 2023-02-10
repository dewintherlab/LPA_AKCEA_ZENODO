#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Function definitions

#-----------------------------------------------------------------------------------------
#Calculate GO and pathway enrichment
#Define a function to plot the GO graphs
plot_GO_from_table <- function(go_data, name="GO_terms.pdf", sig_cut=0.05){
  go_data[,1] <- gsub("_", " ", go_data[,1])
  go_data[,3] <- -log10(go_data[,2])
  colnames(go_data) <- c("term","padj","log10_FDR_qvalue")
  sig_cut <- -log10(sig_cut)
  label_width <- 1 + (max(strwidth(go_data$term, units = "inches")) * 2.54) #Converting to cm
  print(label_width)
  ggplot(go_data,aes(x = reorder(term, log10_FDR_qvalue)  , y = log10_FDR_qvalue)) +
    geom_bar(stat="identity", fill = "Dodgerblue") + 
    coord_flip() +
    ylab("-log10(FDR q-value)") +
    xlab("") +
    #ylim(0,roundCeiling(x = max(go_data$log10_FDR_qvalue))) +
    ylim(0,30) +
    theme(panel.background=element_rect(fill="white"),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm")) +
    geom_hline(yintercept = sig_cut, col = "red", linetype = "dashed")
  ggsave(name,width = (label_width / 2.54) + 7.36)
} #plot_GO_from_table()

#Define helper function for rounding the axis of the plots
roundCeiling <- function(x) {
  if(length(x) != 1) stop("'x' must be of length 1")
  if(x < 5){
    return(5)
  }
  else if(x < 10){
    return(10)
  }
  else if(x < 15){
    return(15)
  }
  else if(x < 100){
    round(x+5,-1)
  }
  else if(x < 1000){
    round(x+50,-2)
  }
  else{
    round(x+500,-3)
  }
} #roundCeiling()

get_ontology <- function(res, organism = "H", name = "treated", return.data = F, padj.cutoff = 0.1, mydir = "differential_expression/", log2FC.cutoff = 0){
  #Get background of all expressed genes in our dataset
  universe <- rpkmData$EntrezGene
  
  #We are using human data
  if (organism == "H"){
    #Up regulated genes
    #Pathways
    #Extract genes from res with under cut off (default 0.1 and 0)
    genes <- data.frame(entrezGene = row.names(subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)), log2FoldChange = subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)$log2FoldChange)
    
    #Make symbols object
    gene_symbols <- data.frame(entrezGene = genes, Symbol = subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)$Symbol)
    
    #Create annotation table
    gs.annots <- buildCustomIdx(genes$entrezGene,
                                species = "human",
                                gsets = ont_sets,
                                label = "CP_and_H",
                                name = "CP_and_H"
    )
    
    #Get pathway enrichments
    gsa.pathways <- egsea.ora(geneIDs = genes$entrezGene,
                              logFC = genes$log2FoldChange,
                              title = "x", 
                              universe = universe,
                              gs.annots = gs.annots,
                              symbolsMap = gene_symbols,
                              sort.by = "p.adj",
                              report.dir = ".",
                              num.threads = 4,
                              report = F
    )
    
    #GO
    #Create annotation table
    gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezGene, species = "human", geneSets = "c5")
    
    #Get GO enrichments
    gsa.go <- egsea.ora(geneIDs = genes$entrezGene,
                        gs.annots = gs.annots,
                        title = "x",
                        universe = universe,
                        symbolsMap = gene_symbols,
                        sort.by = "p.adj",
                        report.dir = ".",
                        num.threads = 4,
                        report = F
    )
    
    #Extract top 20 enriched sets and plot
    top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)), padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
    top.GO <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
    
    plot_GO_from_table(top.pathways, name = paste(mydir, name, ".pathways_UP.pdf", sep = ""))
    plot_GO_from_table(top.GO, name = paste(mydir, name, ".GO_terms_UP.pdf", sep = ""))
    
    #Safe the tables for optional return
    if(return.data){
      up.top.pathways <- top.pathways
      up.top.GO <- top.GO
    }
    
    
    #--------------------------------
    #Down regulated genes
    #Pathways
    #Extract genes from res with padj under cut off (default 0.1)
    genes <- data.frame(entrezGene = row.names(subset(res, padj <padj.cutoff & log2FoldChange < -log2FC.cutoff)), log2FoldChange = subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff)$log2FoldChange)
    
    #Make symbols object
    gene_symbols <- data.frame(entrezGene = genes, Symbol = subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff)$Symbol)
    
    #Create annotation table
    gs.annots <- buildCustomIdx(genes$entrezGene,
                                species = "human",
                                gsets = ont_sets,
                                label = "CP_and_H",
                                name = "CP_and_H"
    )
    
    #Get pathway enrichments
    gsa.pathways <- egsea.ora(geneIDs = genes$entrezGene,
                              logFC = genes$log2FoldChange,
                              title = "x",
                              universe = universe,
                              gs.annots = gs.annots,
                              symbolsMap = gene_symbols,
                              sort.by = "p.adj",
                              report.dir = ".",
                              num.threads = 4,
                              report = F
    )
    
    #GO
    #Create annotation table
    gs.annots <- buildMSigDBIdx(entrezIDs = genes$entrezGene, species = "human", geneSets = "c5")
    
    #Get GO enrichments
    gsa.go <- egsea.ora(geneIDs = genes$entrezGene,
                        gs.annots = gs.annots,
                        title = "x",
                        universe = universe,
                        symbolsMap = gene_symbols,
                        sort.by = "p.adj",
                        report.dir = ".",
                        num.threads = 4,
                        report = F
    )
    
    #Extract top 20 enriched sets and plot
    top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)), padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
    top.GO <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
    
    plot_GO_from_table(top.pathways, name = paste(mydir, name, ".pathways_DOWN.pdf", sep = ""))
    plot_GO_from_table(top.GO, name = paste(mydir, name, ".GO_terms_DOWN.pdf", sep = ""))
    
    #Return the tables if wanted
    if(return.data){
      returnList < c(up.top.pathways,up.top.go,top.pathways,top.GO)
      names(returnList) <- c("Pathways_Up", "GO_terms_Up","Pathways_Down", "GO_terms_Down")
      return(returnList)
    }
    
    #Done with human data
    
    #--------------------------------------------------------------
    
  } else if (organism == "M"){
    #We are using mouse data
    
    #Up regulated genes
    #Pathways
    #Extract genes from res with padj under cut off (default 0.1)
    genes <- subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)$RefSeq_ID
    
    #Convert to human entrezGene IDs
    hum_genes <- getLDS(attributes = "refseq_mrna",
                        filters = "refseq_mrna",
                        values = genes,
                        mart = mouse,
                        attributesL = "entrezgene",
                        martL = human
    )
    hum_universe <- getLDS(attributes = "entrezgene",
                           filters = "entrezgene",
                           values = universe,
                           mart = mouse,
                           attributesL = "entrezgene",
                           martL = human
    )

    #Add Symbols
    symbols <- subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)[,c("RefSeq_ID","Symbol")]
    hum_genes <- merge(hum_genes,symbols,by = "RefSeq_ID")
    
    #Make symbol object
    gene_symbols <- hum_genes[,c("entrezGene","Symbol")]
    
    #Create annotation table
    gs.annots <- buildCustomIdx(hum_genes$entrezGene,
                                species = "human",
                                gsets = ont_sets,
                                label = "CP_and_H",
                                name = "CP_and_H"
    )
    
    #Get pathway enrichments
    gsa.pathways <- egsea.ora(geneIDs = hum_genes$entrezGene,
                              logFC = hum_genes$log2FoldChange,
                              title = "x",
                              universe = hum_universe,
                              gs.annots = gs.annots,
                              symbolsMap = gene_symbols,
                              sort.by = "p.adj",
                              report.dir = file.path(paste(name,"pathways",sep = "_")),
                              num.threads = 4,
                              report = F
    )
    
    #GO
    #Extract genes from res with padj under cut off (default 0.1)
    genes <- row.names(subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff))
    
    #Extract Symbols
    gene_symbols <- cbind(genes,subset(res, padj < padj.cutoff & log2FoldChange > log2FC.cutoff)$Symbol)
    
    #Create annotation table
    gs.annots <- buildMSigDBIdx(entrezIDs = genes, species = "mouse", geneSets = "c5")
    
    #Get GO enrichments
    gsa.go <- egsea.ora(geneIDs = genes,
                        gs.annots = gs.annots,
                        title = "x",
                        universe = hum_universe,
                        symbolsMap = gene_symbols,
                        sort.by = "p.adj",
                        report.dir = file.path(paste(name,"pathways",sep = "_")),
                        num.threads = 4,
                        report = F
    )
    
    #Extract top 20 enriched sets and plot
    top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)), padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
    top.GO <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
    
    plot_GO_from_table(top.pathways, name = paste(name, " pathways UP",".pdf", sep = ""))
    plot_GO_from_table(top.GO, name = paste(name, " GO terms UP", ".pdf", sep = ""))
    
    #Safe the tables for optional return
    if(return.data){
      up.top.pathways <- top.pathways
      up.top.GO <- top.GO
    }
    
    #--------------------------------
    #Down regulated genes
    #Pathways
    #Extract genes from res with padj under cut off (default 0.1)
    genes <- subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff)$RefSeq_ID
    
    #Convert to human entrezGene IDs
    hum_genes <- getLDS(attributes = "refseq_mrna",
                        filters = "refseq_mrna",
                        values = genes,
                        mart = mouse,
                        attributesL = "entrezgene",
                        martL = human
    )
    colnames(hum_genes) <- c("RefSeq_ID","entrezGene")
    
    hum_universe <- getLDS(attributes = "entrezgene",
                           filters = "entrezgene",
                           values = universe,
                           mart = mouse,
                           attributesL = "entrezgene",
                           martL = human
    )
    
    #Add Symbols
    symbols <- subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff)[,c("RefSeq_ID","Symbol")]
    hum_genes <- merge(hum_genes,symbols,by = "RefSeq_ID")
    
    #Make symbol object
    gene_symbols <- hum_genes[,c("entrezGene","Symbol")]
    
    #Create annotation table
    gs.annots <- buildCustomIdx(hum_genes$entrezGene,
                                species = "human",
                                gsets = ont_sets,
                                label = "CP_and_H",
                                name = "CP_and_H"
    )
    
    #Get pathway enrichments
    gsa.pathways <- egsea.ora(geneIDs = hum_genes$entrezGene,
                              logFC = hum_genes$log2FoldChange,
                              title = "x",
                              universe = hum_universe,
                              gs.annots = gs.annots,
                              symbolsMap = gene_symbols,
                              sort.by = "p.adj",
                              report.dir = file.path(paste(name,"pathways",sep = "_")),
                              num.threads = 4,
                              report = F
    )
    
    #GO
    #Extract genes from res with padj under cut off (default 0.1)
    genes <- row.names(subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff))
    
    #Extract Symbols
    gene_symbols <- cbind(genes,subset(res, padj < padj.cutoff & log2FoldChange < -log2FC.cutoff)$Symbol)
    
    #Create annotation table
    gs.annots <- buildMSigDBIdx(entrezIDs = genes, species = "mouse", geneSets = "c5")
    
    #Get GO enrichments
    gsa.go <- egsea.ora(geneIDs = genes,
                        title = "x",
                        universe = hum_universe,
                        gs.annots = gs.annots,
                        symbolsMap = gene_symbols,
                        sort.by = "p.adj",
                        report.dir = file.path(paste(name,"pathways",sep = "_")),
                        num.threads = 4,
                        report = F
    )
    
    #Extract top 20 enriched sets and plot
    top.pathways <- data.frame(name = row.names(topSets(gsa.pathways, number = 20, names.only = F)), padj = topSets(gsa.pathways, number = 20, names.only = F)$p.adj)
    top.GO <- data.frame(name = row.names(topSets(gsa.go, number = 20, names.only = F)), padj = topSets(gsa.go, number = 20, names.only = F)$p.adj)
    
    plot_GO_from_table(top.pathways, name = paste(name, " pathways DOWN",".pdf", sep = ""))
    plot_GO_from_table(top.GO, name = paste(name, " GO terms DOWN", ".pdf", sep = ""))
    
    #Return the tables if wanted
    if(return.data){
      returnList < c(up.top.pathways,up.top.go,top.pathways,top.GO)
      names(returnList) <- c("Pathways_Up", "GO_terms_Up","Pathways_Down", "GO_terms_Down")
      return(returnList)
    }
    
    #Done with mouse data
    
    #-----------------------------------------------------------------------------
    
  } else {
    #Unsupported organism
    
    print("Unrecognised organism!")
    print("Use H for human or M for mouse.")
    
  }
} #get_ontology(res)

#------------------------------------------------------------------------------------------------
#volcano_plot
volcano_plot <- function(res, padj = 0.1, log2FC = 1, outliers = T, labels = F, maxLabels = 10, maxXlim = 5, minXlim = -5, maxYlim = 20, autoScaleAxes = T, autoScaleLabels = T){
  
  #Subset the data on padj and melt for ggplot
  d <- res[, c("Symbol", "pvalue", "padj", "log2FoldChange")]
  d <- d[is.finite(d$log2FoldChange),]
  d <- d[!is.na(d$log2FoldChange),]
  d <- d[is.finite(d$pvalue),]
  d <- d[!is.na(d$pvalue),]
  d <- d[is.finite(d$padj),]
  d <- d[!is.na(d$padj),]
  
  #Set axes limits and remove outliers if wanted
  if(autoScaleAxes){
    if(outliers){
      maxYlim <- roundCeiling(max(-log10(d$pvalue)))
      
      maxXlim <-  max(d$log2FoldChange)
      minXlim <- -roundCeiling(-min(d$log2FoldChange))
      maxXlim <- max(c(maxXlim,-minXlim))
      minXlim <- -maxXlim
    }else{
      maxYlim <- roundCeiling(quantile(-log10(d$pvalue), prob = 0.99))
      
      maxXlim <- roundCeiling(max(d[-log10(d$pvalue) < maxYlim,]$log2FoldChange))
      minXlim <- -roundCeiling(-min(d[-log10(d$pvalue) < maxYlim,]$log2FoldChange))
      maxXlim <- max(c(maxXlim,-minXlim))
      minXlim <- -maxXlim
      
      d <- d[which(-log10(d$pvalue) < maxYlim),]
    }
  }
  
  #Set color coding
  d$color <- "black"
  d[which(d$padj < padj), "color"] <- "orange"
  d[which(abs(d$log2FoldChange) > log2FC), "color"] <- "red"
  d[which(d$padj < padj & abs(d$log2FoldChange) > log2FC), "color"] <- "green"
  d$color <- factor(d$color, levels = c("black", "green", "orange", "red"))
  
  #set vars for use in gg call
  padj_Var <- padj
  thelim <- list()
  
  #and plot!
  if(labels){
    #Get padj and FC cut-off for labels
    if(autoScaleLabels){
      padj.FC.lims <- get.padj.FC.lims(res)
      thelim       <- get.opt.padj.FC.lim(maxLabels,padj.FC.lims) #Plot only < n observations
    }else{
      thelim[1] <- log2FC
      thelim[2] <- padj
    }
    
    ggplot(d) +
      geom_point(aes(x = log2FoldChange, y = -log10(padj), color = color)) +
      scale_color_manual(name = "Legend", 
                         values = setNames(
                           c("slategray3", "springgreen", "goldenrod", "firebrick"), 
                           c("black", "green", "orange", "red")),
                         labels = c("ns", 
                                    substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC, " & Adjusted P-value < ", padj_Var), list(log2FC=log2FC, padj_Var=padj_Var)), 
                                    paste("Adjusted P-value <", padj_Var, sep = " "), 
                                    substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC), list(log2FC=log2FC))),
                         drop = F) +
      geom_label_repel(
        aes(x = log2FoldChange, y = -log10(padj), label = ifelse(padj < thelim[2] & abs(log2FoldChange) > thelim[1], as.character(Symbol),'')),
        box.padding = 0.35, point.padding = 0.5,
        segment.color = 'grey50') +
      geom_vline(xintercept = 0, size = 0.05, color = "red") +
      geom_vline(xintercept = log2FC, size = 0.05, color = "red") +
      geom_vline(xintercept = -log2FC, size = 0.05, color = "red") +
      geom_hline(yintercept = -log10(padj), size = 0.05, color = "red") +
      xlim(c(minXlim, maxXlim)) +
      xlab(expression(log[2] ~("Fold Change"))) +
      ylab(expression(-log[10] ~("adjusted P-value"))) +
      scale_y_continuous(expand = c(0, 0.1), limits = c(0, maxYlim) ) +
      theme_light() +
      theme(axis.text = element_text(size = 16),
            axis.ticks = element_line(size = 0.5),
            axis.title = element_text(size = 16),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_line(size = 1),
            panel.border = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 14)
            
      )
  }else{
    ggplot(d) +
      geom_point(aes(x = log2FoldChange, y = -log10(padj), color = color)) +
      scale_color_manual(name = "Legend", 
                         values = setNames(
                           c("slategray3", "springgreen", "goldenrod", "firebrick"), 
                           c("black", "green", "orange", "red")),
                         labels = c("ns", 
                                    substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC, " & Adjusted P-value < ", padj_Var), list(log2FC=log2FC, padj_Var=padj_Var)), 
                                    paste("Adjusted P-value <", padj_Var, sep = " "), 
                                    substitute(paste(group("|", log[2]("FC"), "|"), " > ", log2FC), list(log2FC=log2FC))),
                         drop = F) +
      geom_vline(xintercept = 0, size = 0.05, color = "red") +
      geom_vline(xintercept = log2FC, size = 0.05, color = "red") +
      geom_vline(xintercept = -log2FC, size = 0.05, color = "red") +
      geom_hline(yintercept = -log10(padj), size = 0.05, color = "red") +
      xlim(c(minXlim, maxXlim)) +
      xlab(expression(log[2] ~("Fold Change"))) +
      ylab(expression(-log[10] ~("adjusted P-value"))) +
      scale_y_continuous(expand = c(0, 0.1), limits = c(0, maxYlim) ) +
      theme_light() +
      theme(axis.text = element_text(size = 16),
            axis.ticks = element_line(size = 0.5),
            axis.title = element_text(size = 16),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_line(size = 1),
            panel.border = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 14)
      )
  }
}

#-----------------------------------------------------------------------------------------
#Find number of observations for a padj and FC combination
num.padj.FC <- function(p,FC,res=res){
  length(which(res$padj <= p & (res$log2FoldChange > FC | res$log2FoldChange < -FC)))
} #num.padj.FC()

#-----------------------------------------------------------------------------------------
#Create matrix of yields of padj and FC combinations
get.padj.FC.lims <- function(res=res){
  ps <- c(0.1,0.05,0.01,0.001,0.0001,0.00001)
  FCs <- c(0.5,1,2,3,4,5,6,7,8,9,10)
  x <- data.frame()
  
  for(i in 1:length(ps)){
    for(j in 1:length(FCs)){
      x[i,j] <- num.padj.FC(ps[i],FCs[j],res)
    }
  }
  colnames(x) <- FCs
  row.names(x) <- ps
  return(x)
} #get.padj.FC.lims()


#-----------------------------------------------------------------------------------------
#Find the optimum (lowest FC) padj and FC combination that stays below n observations
get.opt.padj.FC.lim <- function(n,m){
  for (i in 1:length(colnames(m))){
    x <- which(m[,i] < n)
    if (length(x) > 0){
      return(c(as.numeric(colnames(m)[i]),as.numeric(row.names(m)[min(x)])))
      break
    }
  }
} #get.opt.padj.FC.lim()


#-----------------------------------------------------------------------------------------
#View box plots of single genes within R
ggPlotCounts <- function(theGene,intgroup="treatment",subgroup=1:length(colnames(dds)), res.obj = res, dds.obj = dds, colData.obj = colData){
  print(theGene)
  print(res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"])
  d <- plotCounts(dds.obj[,subgroup], gene=res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"], intgroup=intgroup,returnData=T)
  ggplot(d, aes_string(intgroup, "count")) +
    geom_point(position=position_jitter(w=0.1, h=0),mapping = aes(color=colData.obj[subgroup,intgroup])) +
    scale_colour_manual(name=intgroup, values = ann_colors$Group, breaks = names(ann_colors$Group), labels = c("Healthy controls - Baseline LPA level",
                                                                             "Healthy controls - High LPA level",
                                                                             "Patients - High LPA level",
                                                                             "Patients - Treatment reduced LPA level")) +
    ggtitle(theGene) +
    ylab("Normalized count") +
    xlab("") +
    stat_summary(fun.y = median,
                 fun.ymin = function(z) {quantile(z, 0.25)},
                 fun.ymax = function(z) {quantile(z, 0.75)},
                 geom = "crossbar",
                 width = 0.5,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z,0.25)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
                 fun.ymax = function(z) {quantile(z,0.95)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z, 0.05)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
                 fun.ymax = function(z) {quantile(z, 0.95)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,intgroup])
    ) +
    theme_light() +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          aspect.ratio = 2/1, 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "right"
    )
} #ggPlotCounts()

ggPlotCountsGR <- function(theGene,intgroup="treatment",subgroup=1:length(colnames(dds)), res.obj = res, dds.obj = dds, colData.obj = colData){
  print(theGene)
  print(res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"])
  d <- plotCounts(dds.obj[,subgroup], gene=res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"], intgroup=intgroup,returnData=T)
  ggplot(d, aes_string(intgroup, "count")) +
    geom_point(position=position_jitter(w=0.1, h=0), size = 3, aes(color = colData.obj$GR)) +
    scale_colour_manual(name="GR", values = ann_colors$GR, labels = c("50-90% reduction - before treatment",
                                                                      "10-50% reduction - before treatment", 
                                                                      "50-90% reduction - after treatment",
                                                                      "10-50% reduction - after treatment")) +
    ggtitle(theGene) +
    ylab("Normalized count") +
    xlab("") +
    stat_summary(fun.y = median,
                 fun.ymin = function(z) {quantile(z, 0.25)},
                 fun.ymax = function(z) {quantile(z, 0.75)},
                 geom = "crossbar",
                 width = 0.5,
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z,0.25)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
                 fun.ymax = function(z) {quantile(z,0.95)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z, 0.05)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
                 fun.ymax = function(z) {quantile(z, 0.95)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    theme_light() +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          aspect.ratio = 2/1, 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "right"
    )
}

ggPlotCountsRes <- function(theGene,intgroup="treatment",subgroup=1:length(colnames(dds)), res.obj = res, dds.obj = dds, colData.obj = colData){
  print(theGene)
  print(res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"])
  d <- plotCounts(dds.obj[,subgroup], gene=res.obj[which(res.obj$Symbol == theGene),"RefSeq_ID"], intgroup=intgroup,returnData=T)
  ggplot(d, aes_string(intgroup, "count")) +
    geom_point(position=position_jitter(w=0.1, h=0), size = 3, aes(color = colData.obj$GR)) +
    scale_colour_manual(name="GR", values = ann_colors$GR, labels = c("50-90% reduction - after treatment",
                                                                      "10-50% reduction - after treatment")) +
    ggtitle(theGene) +
    ylab("Normalized count") +
    xlab("") +
    stat_summary(fun.y = median,
                 fun.ymin = function(z) {quantile(z, 0.25)},
                 fun.ymax = function(z) {quantile(z, 0.75)},
                 geom = "crossbar",
                 width = 0.5,
                 size = 0.5,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z,0.25)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
                 fun.ymax = function(z) {quantile(z,0.95)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z, 0.05)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
                 fun.ymax = function(z) {quantile(z, 0.95)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color = colData(dds.obj)[subgroup,"GR"])
    ) +
    theme_light() +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          aspect.ratio = 2/1, 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "right"
    )
}

#-----------------------------------------------------------------------------------------
#View box plots of single genes within R
ggPlotCountsInteractions <- function(theGene,intgroup=c("time","condition")){
  print(theGene)
  print(res[which(res$Symbol == theGene),"RefSeq_ID"])
  d <- plotCounts(dds[,], gene=res[which(res$Symbol == theGene),"RefSeq_ID"], intgroup=intgroup,returnData=T)
  ggplot(d, aes(x = get(intgroup[1]), y = count)) +
    geom_point(position=position_jitter(w=0.1, h=0),mapping = aes(color= inhibition)) +
    scale_colour_manual(name=intgroup[2], values = c("dodgerblue", "coral")) +
    ggtitle(theGene) +
    ylab("Normalized count") +
    xlab("") +
    stat_summary(fun.y = median,
                 fun.ymin = function(z) {quantile(z, 0.25)},
                 fun.ymax = function(z) {quantile(z, 0.75)},
                 geom = "crossbar",
                 width = 0.5,
                 size = 0.1,
                 mapping = aes(color = get(intgroup[2]))
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z,0.25)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color =   get(intgroup[2]))
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.75)},
                 fun.ymax = function(z) {quantile(z,0.95)},
                 geom = "linerange",
                 size = 0.1,
                 mapping = aes(color =   get(intgroup[2]))
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.05)},
                 fun.ymax = function(z) {quantile(z, 0.05)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color =   get(intgroup[2]))
    ) +
    stat_summary(fun.ymin = function(z) {quantile(z, 0.95)},
                 fun.ymax = function(z) {quantile(z, 0.95)},
                 geom = "errorbar",
                 width = 0.25,
                 size = 0.1,
                 mapping = aes(color =   get(intgroup[2]))
    ) +
    theme_light() +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          legend.position = "right"
    )
} #ggPlotCountsInteractions()

#-----------------------------------------------------------------------------------------
#PCA Plot
PCA.plot <- function(rld,intgroup="Sex"){
  mat.pca <- as.data.frame(prcomp(t(assay(rld)))$x)
  mat.pca <- cbind(mat.pca, intgroup=colData[,intgroup]) #Add metadata
  percentVar <- round(100 * prcomp(t(assay(rld)))$sdev^2/sum(prcomp(t(assay(rld)))$sdev^2))
  
  ggplot(mat.pca, aes(PC1, PC2, color = intgroup)) +
    geom_point(size = 2) +
    scale_colour_manual(name=intgroup, values = c("dodgerblue", "coral")) +
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
}
