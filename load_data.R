## install.packages("BiocManager")
## BiocManager::install(c("DESeq2", "sva", "RColorBrewer", "gplots", "BiocParallel", "genefilter","ggplot2","grid", "scatterplot3d", "pheatmap", "reshape2", "org.Hs.eg.db", "edgeR", "EGSEA", "ggrepel", "biomaRt", "rafalib"), ask = F)

#Load packages
library(DESeq2)
library(sva)
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
library(ggrepel)
library(biomaRt)
library(rafalib)
register(MulticoreParam(4))
source("functions.R")


#-----------------------------------------------------------------------------------------
#Set up gene ontology objects for pathway analyses
ont_sets <- readList("cp_and_h.gmt.txt")

marts <- listMarts()
mart  <- marts[1,1]

human <- useMart(mart, dataset = "hsapiens_gene_ensembl")

#------------------------------------------------------------------------------------------
#Load Data
#AKCEA
countData <- read.table("raw_read.quant_table.clean.txt",header=T,row.names=1)
rpkmData  <- read.table("RPKM.quant_table.clean.txt",header=T,row.names=1)
colData   <- read.table("metadata.table.txt",header=T,row.names=1, colClasses = "factor")

colData$Patient   <- relevel(colData$Patient, "2")
colData$Visit     <- relevel(colData$Visit, "1")
colData$Treatment <- relevel(colData$Treatment, "placebo")
colData$Group     <- relevel(colData$Group, "baseline")
colData$Response     <- relevel(colData$Response, "Low")

#Remove outlier S08
#rpkmData  <- rpkmData[,which(colnames(rpkmData)   != "S08")]
#countData <- countData[,which(colnames(countData) != "S08")]
#colData   <- colData[which(row.names(colData)     != "S08"),]

#Remove placebo samples and matching baseline controls
countData <- countData[,grep("placebo",colData$Group, invert = T)]
rpkmData  <- rpkmData[,grep("placebo",colData$Group, invert = T)]
colData   <- colData[grep("placebo",colData$Group, invert = T),]

#Baseline LPA
#Load Data
LPA.countData <- read.table("LPA.raw_read_quant.table.clean.txt",header=T,row.names=1)
LPA.rpkmData  <- read.table("LPA.RPKM_quant.table.clean.txt",header=T,row.names=1)
LPA.colData   <- read.table("LPA.metadata.table.txt",header=T,row.names=1, colClasses = "factor")

LPA.colData$LPA  <- relevel(LPA.colData$LPA, "Low")
LPA.colData$Sex  <- relevel(LPA.colData$Sex, "F")

#Remove outliers S01, S02, S03, S09, S17, S11
LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S01")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S01")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S01"),]

LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S02")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S02")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S02"),]

LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S03")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S03")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S03"),]

LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S09")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S09")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S09"),]

LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S17")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S17")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S17"),]

LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S11")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S11")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S11"),]

#Remove erronous samples P63 and P67 (S02, S06)
#S02 was already removed
LPA.rpkmData  <- LPA.rpkmData[,which(colnames(LPA.rpkmData)   != "S06")]
LPA.countData <- LPA.countData[,which(colnames(LPA.countData) != "S06")]
LPA.colData   <- LPA.colData[which(row.names(LPA.colData)     != "S06"),]

#Rectify wrong classification of P73 (S11)
#Already removed as outlier
#colData["S11","LPA"] <- "High"


#-----------------------------------------------------------------------------------------
#Merge data sets
#Remove superfluous metadata
colData$Visit     <- NULL
colData$Treatment <- NULL

LPA.colData$Group     <- LPA.colData$LPA
LPA.colData$Response  <- "none"
LPA.colData$LPA       <- NULL
LPA.colData$Sex       <- NULL
LPA.colData$Patient   <- "NA"

#Adjust row and col names so they don't clash
colnames(countData) <- paste("Patient",colnames(countData),sep = "_")
colnames(rpkmData)  <- paste("Patient",colnames(rpkmData),sep = "_")
row.names(colData)  <- paste("Patient",row.names(colData),sep = "_")

colnames(LPA.countData) <- paste("Healthy",colnames(LPA.countData),sep = "_")
colnames(LPA.rpkmData)  <- paste("Healthy",colnames(LPA.rpkmData),sep = "_")
row.names(LPA.colData)  <- paste("Healthy",row.names(LPA.colData),sep = "_")
colnames(LPA.colData)   <- c("Flowcell", "Group", "Response", "Patient")

#Merge the tables
colData <- rbind(colData,LPA.colData)

countData            <- merge(countData,LPA.countData, by = 0)
row.names(countData) <- countData$Row.names
countData$Row.names  <- NULL

rpkmData            <- merge(rpkmData,LPA.rpkmData, by = 0)
row.names(rpkmData) <- rpkmData$Row.names
rpkmData$Row.names  <- NULL

#Rename the groups
colData$Group <- as.character(colData$Group)

colData[which(colData$Group == "Low"      ), "Group"] <- "Healthy_baseline_LPA"
colData[which(colData$Group == "High"     ), "Group"] <- "Healthy_high_LPA"
colData[which(colData$Group == "baseline" ), "Group"] <- "Patient_high_LPA"
colData[which(colData$Group == "treatment"), "Group"] <- "Patient_treated_LPA"

colData$Group <- as.factor(colData$Group)
colData$Group <- relevel(colData$Group, "Healthy_baseline_LPA")

colData$Flowcell <- as.factor(colData$Flowcell)
colData$Patient <- as.factor(colData$Patient)

#Define column numbers
HEALTHY_BASELINE <- which(colData$Group == "Healthy_baseline_LPA")
HEALTHY_HIGH     <- which(colData$Group == "Healthy_high_LPA")
PATIENT_HIGH     <- which(colData$Group == "Patient_high_LPA")
PATIENT_TREATED  <- which(colData$Group == "Patient_treated_LPA")
META_COLS        <- (length(colnames(rpkmData))+1):(length(colnames(rpkmData))+3) #WE will later add 3 columns to the rpkm table


#-----------------------------------------------------------------------------------------
#Remove low or high count genes
#Keep median RPKM > 0 in all samples
keep      <- apply(rpkmData,1,sum) > 0
rpkmData  <- rpkmData[keep,]
countData <- countData[keep,]

#Remove median RPKM > 2000
keep      <- apply(rpkmData,1,median) < 2000
rpkmData  <- rpkmData[keep,]
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
rpkmData  <- rpkmData[o,]
d <- duplicated(rpkmData$Symbol)
countData <- countData[!d,]
rpkmData  <- rpkmData[!d,]


#-----------------------------------------------------------------------------------------
#Remove ribosomal, mitochondrial, non coding, and MIR genes
#keep <- grep("^MT-|^RPL|^RPS|^RNU|^RNVU|^MIR|^SNOR|^SCAR|^LINC|^LOC", rpkmData$Symbol, invert = T)
#rpkmData <- rpkmData[keep,]
#countData <- countData[keep,]

#-----------------------------------------------------------------------------------------
#Deduce Sex
colData$Sex <- rpkmData[which(rpkmData[,"Symbol"] == "UTY"), 1:(length(colnames(rpkmData))-3), drop = T] == 0
colData[which(colData$Sex == FALSE ), "Sex"] <- "Male"
colData[which(colData$Sex == TRUE), "Sex"] <- "Female"
colData$Sex <- as.factor(colData$Sex)
colData$Sex <- relevel(colData$Sex, "Female")

#------------------------------------------------------------------------------------------
#Create blocking group
colData$Blocks <- as.factor(paste(colData$Flowcell, colData$Sex, sep = "."))

#------------------------------------------------------------------------------------------
#Setup heatmap annotation colors
ann_colors <- list(
  Group = c(Healthy_baseline_LPA = "chartreuse3", Healthy_high_LPA = "coral", Patient_high_LPA = "firebrick1", Patient_treated_LPA = "chartreuse4"),
  Flowcell = c(HCLLKDSXX = "dodgerblue", HCMKVDSXX = "dodgerblue4"),
  Sex = c(Male = "cornflowerblue", Female = "darksalmon"),
  Response = c(Low = "maroon4", High = "maroon1"),
  GR = c(Patient_high_LPA.High = "firebrick1", Patient_high_LPA.Low = "firebrick4", Patient_treated_LPA.High = "green", Patient_treated_LPA.Low = "green4")
)


#Setup the data object fo rthe various comparisons
source("setup_data_all_groups.R")
source("setup_data_AvB.R")
source("setup_data_AvC.R")
source("setup_data_CvD.R")
source("setup_data_stratified.CvD.R")
source("setup_data_ultra_stratified.CvD.R")


#Run the DE analyses
source("DE Code.R")
source("AvB.DE_Code.R")
source("AvC.DE_Code.R")
source("CvD.DE_Code.R")
source("stratified.CvD.DE_Code.R")
source("ultra_stratified.CvD.DE_Code.R")

#Prepare Data fro GSEA analysis
GSEAData <- rpkmData[which(!is.na(rpkmData$Symbol)),1:(length(colnames(rpkmData))-3)]
row.names(GSEAData) <- rpkmData[which(!is.na(rpkmData$Symbol)),"Symbol"]
write.table(transform(GSEAData, Symbol = rownames(GSEAData))[,c(ncol(GSEAData)+1,1:length(colnames(GSEAData)))],file="GSEA_data.txt", row.names = F, quote=F, sep="\t")

write.table(colData[,"Group"],file = "GSEA_phenotype_data.cls", quote = F, row.names = F, col.names = F, sep = " ")
