########################
### Install packages ###
########################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("TCGAbiolinks","SummarizedExperiment","edgeR","org.Hs.eg.db","GSEABase", "ExpressionNormalizationWorkflow","GSVA","monocle","ggplot2"))

################################
### Load TCGA-LIHC-FPKM data ###
################################

## Set your working directory (or use Control + Shift + H to choose it) ##
## The working directory should contain the RData files "181208-TCGA-LIHC-hg38-FPKM.rda" and "TCGA-Monocle.RData" ##
setwd("...")


### prepare data into SummarizedExperiment

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

load("181208-TCGA-LIHC-hg38-FPKM.rda")
load("TCGA-Monocle.RData")

LIHC_FPKM <- data
rm(data)
dim(LIHC_FPKM)
#View(LIHC_FPKM)
#View(assay(LIHC_FPKM))
#View(colData(LIHC_FPKM))


### Creat gene_annotation
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
gene_annotation <- as.data.frame(rownames(assay(LIHC_FPKM)))
rownames(gene_annotation) <- gene_annotation[,1]
colnames(gene_annotation)[1] <- "ENSEMBL"
gene_annotation$gene_short_name <- mapIds(org.Hs.eg.db,
                          keys=as.character(rownames(gene_annotation)),
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
gene_annotation$ENTREZID <- mapIds(org.Hs.eg.db,
                            keys=as.character(rownames(gene_annotation)),
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
genes_have_sym <- rownames(gene_annotation[!is.na(gene_annotation$gene_short_name),])
length(genes_have_sym)
#View(gene_annotation)


##########################################
### Construct MRS associated gene sets ###
##########################################
library(GSEABase)
#OutputDEG <- read.csv("Differential expression analysis result.csv",
#                      header=TRUE, as.is=TRUE)

#View(OutputDEG)
upDEG <- OutputDEG[OutputDEG$adj.P.Val<=0.01 & OutputDEG$logFC>1,]
upDEG$EntrezID <- as.character(upDEG$EntrezID)
downDEG <- OutputDEG[OutputDEG$adj.P.Val<=0.01 & OutputDEG$logFC<(-1),]
downDEG$EntrezID <- as.character(downDEG$EntrezID)
dim(upDEG)
dim(downDEG)
mrsUpSet <- GeneSet(unique(upDEG$EntrezID),
                    setName="MRSup",
                    geneIdType=EntrezIdentifier("org.Hs.eg.db"))

mrsDnSet <- GeneSet(unique(downDEG$EntrezID),
                    setName="MRSdn",
                    geneIdType=EntrezIdentifier("org.Hs.eg.db"))

######################################################
### Construct a CD8 exhaustion associated gene set ###
######################################################
#Exhaust <- read.csv("Zheng et al (2017)-T cell exhaustion-related genes.csv",
#                     header=TRUE, as.is=TRUE)
#View(Exhaust)
#dim(Exhaust)
Exhaust$ont <- factor(Exhaust$ont)
Exhaust$gene <- as.character(Exhaust$gene)

exGS <- GeneSet(unique(Exhaust$gene),
                setName="Exhaustion_related_genes",
                geneIdType=EntrezIdentifier("org.Hs.eg.db"))
ImmSets<- GeneSetCollection(c(mrsUpSet,mrsDnSet,exGS)) 
toGmt(ImmSets, "ImmSets.gmt")




#################################
### Construct a ExpressionSet ###
#################################
library(ExpressionNormalizationWorkflow)
#TCGApheno <-  read.delim("181208-TCGA_LIHC_FPKM_Clinical.txt",
#                         header=TRUE, as.is=TRUE,sep = "\t",
#                         na.strings = "NA")
rownames(TCGApheno) <- TCGApheno$barcode


### transform Ensembl IDs into EntrezIDs
keytypes(org.Hs.eg.db)
hccData <- as.data.frame(assay(LIHC_FPKM))
hccData$ENTREZID <- mapIds(org.Hs.eg.db,
                         keys=as.character(rownames(hccData)),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")

## exclude genes without Gene Symbols
hccData <- hccData[!is.na(hccData$ENTREZID),]
hccDataSym <- aggregate(x = hccData[,1:(ncol(hccData)-1)], 
                        by = list(hccData$ENTREZID), FUN = max)
#View(hccDataSym)
rownames(hccDataSym) <- hccDataSym[,1]
hccDataSym[,1] <- NULL

### Construct a ExpressionSet ###
all(colnames(as.matrix(hccDataSym))==rownames(TCGApheno))
TCGAhcc <- expSetobj(as.matrix(hccDataSym), TCGApheno)
annotation(TCGAhcc) <- "org.Hs.eg.db"

##########################
### Single sample GSEA ###
##########################
library(GSVA)
### load MRS-related gene set collection ###
ImmSets <- getGmt("ImmSets.gmt")
ImmSets <- mapIdentifiers(ImmSets, EntrezIdentifier())


### GSVA ###

Imm_ssGSEA <- exprs(gsva(TCGAhcc, ImmSets, min.sz=1, max.sz=Inf,method="gsva",
                         kcdf="Poisson", mx.diff=TRUE))
rownames(Imm_ssGSEA)[1:3] <- paste0(rownames(Imm_ssGSEA)[1:3],"_gsva")
View(Imm_ssGSEA)


## merge pheno data
Imm_ssGSEA <- as.data.frame(t(Imm_ssGSEA))
rownames(Imm_ssGSEA) <- gsub("\\.","-",rownames(Imm_ssGSEA))
View(Imm_ssGSEA)
Imm_ssGSEA$barcode <- rownames(Imm_ssGSEA)
#View(Imm_ssGSEA)
LIHC_ImmSSGSEA <- merge(TCGApheno,Imm_ssGSEA, by="barcode",sort = FALSE,
                        all.x=TRUE)
rownames(LIHC_ImmSSGSEA) <- LIHC_ImmSSGSEA$barcode
View(LIHC_ImmSSGSEA)


################################################
### load TIMER estimation results (optional) ###
################################################

### load TIMER estimated immune infiltrates
#TCGA_TIMER <- as.data.frame(read.csv("TCGA-TIMER.csv",
#                                     header=TRUE, as.is=TRUE))

TCGA_TIMER$SAMPLE <- TCGA_TIMER$barcode
#View(TCGA_TIMER)

LIHC_ImmSSGSEA$SAMPLE <- substr(LIHC_ImmSSGSEA$barcode,1,15)
#View(LIHC_ImmSSGSEA)

### transform TIMER "barcode" (sample code) into real "barcode"
table(LIHC_ImmSSGSEA$SAMPLE %in% TCGA_TIMER$SAMPLE)
LIHC_TIMER <- merge(LIHC_ImmSSGSEA,TCGA_TIMER, by="SAMPLE", sort = FALSE,
                    all.y=FALSE)
colnames(LIHC_TIMER)[1] <- "barcode"
rownames(LIHC_TIMER) <- LIHC_TIMER$barcode.x
LIHC_TIMER[,1] <- LIHC_TIMER$barcode.x
LIHC_TIMER$barcode.y <- NULL
LIHC_TIMER$barcode.x <- NULL
dim(LIHC_TIMER)
#View(LIHC_TIMER)

write.table(as.data.frame(LIHC_TIMER), 
            file="TCGA_LIHC_FPKM_Clinical_cpx.txt",
            quote = FALSE, sep = "\t",row.names=F)


#########################################################
### Convert the SummarizedExperiment into CellDataSet ###
#########################################################
library(monocle)
sample_sheet <-  as.data.frame(read.delim("TCGA_LIHC_FPKM_Clinical_cpx.txt",
                                          header=TRUE, as.is=TRUE,sep = "\t",
                                          na.strings = "NA"))
rownames(sample_sheet) <- sample_sheet[,1]
sample_sheet <- sample_sheet[order(rownames(sample_sheet)),]
#View(sample_sheet)
expr_matrix <- as.matrix(assay(LIHC_FPKM))
expr_matrix <- expr_matrix[,order(colnames(expr_matrix))]
all(colnames(expr_matrix)==rownames(sample_sheet))

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
### Choosing a distribution for your data
cdsLIHC <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily=tobit())



#####################################
### Filtering low-quality samples ###
#####################################
cdsLIHCfilt <- detectGenes(cdsLIHC, min_expr = 0.1)

### filter genes with Symbols and expressed in over 50 samples
expressed_genes <- genes_have_sym[genes_have_sym %in%
                                  row.names(subset(fData(cdsLIHCfilt),
                                    num_cells_expressed >= 50))]
length(expressed_genes)

#############################################
### Calculate ralative cytotoxic activity ###
#############################################

### define a fucntion to calculate geometric mean
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

### find GZMA, PRF1, CD8A, CD8B ids 
GZMA_id <- row.names(subset(fData(cdsLIHCfilt), gene_short_name == "GZMA"))
GZMB_id <- row.names(subset(fData(cdsLIHCfilt), gene_short_name == "GZMB"))
PRF1_id <- row.names(subset(fData(cdsLIHCfilt), gene_short_name == "PRF1"))
CD8A_id <- row.names(subset(fData(cdsLIHCfilt), gene_short_name == "CD8A"))
CD8B_id <- row.names(subset(fData(cdsLIHCfilt), gene_short_name == "CD8B"))

### calculate relative GZMA, PRF1 and CYT expression
pData(cdsLIHCfilt)$GZMA <- exprs(cdsLIHCfilt)[GZMA_id,]
pData(cdsLIHCfilt)$GZMB <- exprs(cdsLIHCfilt)[GZMB_id,]
pData(cdsLIHCfilt)$PRF1 <- exprs(cdsLIHCfilt)[PRF1_id,]
pData(cdsLIHCfilt)$CD8A <- exprs(cdsLIHCfilt)[CD8A_id,]
pData(cdsLIHCfilt)$CD8B <- exprs(cdsLIHCfilt)[CD8B_id,]
pData(cdsLIHCfilt)$CYT <- apply(pData(cdsLIHCfilt)[,c("GZMA","PRF1")],1,gm_mean)
pData(cdsLIHCfilt)$CYT2 <- apply(pData(cdsLIHCfilt)[,c("GZMA","GZMB","PRF1")],1,gm_mean)

pData(cdsLIHCfilt)$relGZMA <- pData(cdsLIHCfilt)$GZMA/apply(pData(cdsLIHCfilt)[,c("CD8A","CD8B")],1,gm_mean)
pData(cdsLIHCfilt)$relPRF1 <- pData(cdsLIHCfilt)$PRF1/apply(pData(cdsLIHCfilt)[,c("CD8A","CD8B")],1,gm_mean)
pData(cdsLIHCfilt)$relCYT <- apply(pData(cdsLIHCfilt)[,c("GZMA","PRF1")],1,gm_mean)/apply(pData(cdsLIHCfilt)[,c("CD8A","CD8B")],1,gm_mean)
pData(cdsLIHCfilt)$relCYT2 <- apply(pData(cdsLIHCfilt)[,c("GZMA","GZMB","PRF1")],1,gm_mean)/apply(pData(cdsLIHCfilt)[,c("CD8A","CD8B")],1,gm_mean)
pData(cdsLIHCfilt)$relExh <- pData(cdsLIHCfilt)$Exhaustion_related_genes_gsva/apply(pData(cdsLIHCfilt)[,c("CD8A","CD8B")],1,gm_mean)
#View(pData(cdsLIHCfilt))



##########################################
### Classifying samples by sample type ###
##########################################
#View(pData(cdsLIHCfilt))
pData(cdsLIHCfilt)$SampleType <- factor(pData(cdsLIHCfilt)$definition) # read sample type info
table(pData(cdsLIHCfilt)$SampleType)

pie <- ggplot(pData(cdsLIHCfilt),
              aes(x = factor(1), fill = SampleType)) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#####################################
### Subset the CellDataSet object ###
#####################################
cdsLIHC_HCC <- cdsLIHCfilt[,pData(cdsLIHCfilt)$SampleType=="Primary solid Tumor"]
table(pData(cdsLIHC_HCC)$SampleType)


################################################
### Constructing MRS associated Trajectories ###
################################################

#OutputDEG <- read.csv("Differential expression analysis result.csv",
#                      header=TRUE, as.is=TRUE)

#View(OutputDEG)
upDEG <- OutputDEG[OutputDEG$adj.P.Val<=0.01 & OutputDEG$logFC>1,]
upGenes <- as.character(upDEG$Symbol)
downDEG <- OutputDEG[OutputDEG$adj.P.Val<=0.01 & OutputDEG$logFC<(-1),]
dnGenes <- as.character(downDEG$Symbol)

mrsGenes <- c(upGenes,dnGenes)
table(duplicated(mrsGenes))
MRS_id <- row.names(subset(fData(cdsLIHC_HCC), gene_short_name %in% mrsGenes))
length(MRS_id)

#ordering_genes <- MRS_id[MRS_id %in% expressed_genes]
ordering_genes <- MRS_id

length(ordering_genes)

cdsLIHC_HCC <- setOrderingFilter(cdsLIHC_HCC, ordering_genes)
cdsLIHC_HCC <- reduceDimension(cdsLIHC_HCC, max_components = 2,
                            method = 'DDRTree')
cdsLIHC_HCC <- orderCells(cdsLIHC_HCC,reverse = FALSE)

colorBrew1 <- c("#43a2ca","#66c2a5","#b3de69","#fdb462","#fb8072")


plot_cell_trajectory(cdsLIHC_HCC, color_by = "State",show_branch_points=F,
                     theta=180,cell_size=2) +
  scale_colour_manual(values=colorBrew1)
## save the plot
ggsave("Plots/MRS_trajectory.pdf", width=100,height=30,units="mm")
#View(pData(cdsLIHC_HCC))

write.table(as.data.frame(pData(cdsLIHC_HCC)), 
            file="TCGA_LIHC_Monocle.txt",
            quote = FALSE, sep = "\t",row.names=F)


  
ggplot(pData(cdsLIHC_HCC), aes(x=State,y=MRSup_gsva),colours=State) +
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=1,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/MRSup.pdf", width=6.5,height=4.5,units="cm")

ggplot(pData(cdsLIHC_HCC), aes(x=State,y=MRSdn_gsva),colours=State) +
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=1,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/MRSdn.pdf", width=6.5,height=4.5,units="cm")

ggplot(pData(cdsLIHC_HCC), aes(x=State,y=CD8_Tcell),colours=State) +
  ylim(0.1,0.3)+
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=2,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/CD8Tcell.pdf", width=6.5,height=4.5,units="cm")


ggplot(pData(cdsLIHC_HCC), aes(x=State,y=Neutrophil),colours=State) +
  ylim(0.05,0.2)+
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=2,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/Neutrophil.pdf", width=6.5,height=4.5,units="cm")

ggplot(pData(cdsLIHC_HCC), aes(x=State,y=relCYT),colours=State) +
  ylim(0,15)+
  scale_y_continuous(breaks = c(0,5,10,15),labels = c("0.00","5","10","15"),
                     limits = c(0,15))+
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=1,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/CYT.pdf", width=6.5,height=4.5,units="cm")



ggplot(pData(cdsLIHC_HCC), aes(x=State,y=Exhaustion_related_genes_gsva),colours=State) +
  #ylim(-0.25,1)+
  scale_colour_manual(values=colorBrew1) +
  geom_violin(aes(colour=State),trim=F,adjust=1,fill="white",
              size=.5)+
  geom_boxplot(width=.15,fill=NA,size=.5,notch=T,
               outlier.colour =NA,colour=colorBrew1)+
  stat_summary(fun.y=mean,geom="point",fill=colorBrew1,
               shape=21,size=.75,colour=colorBrew1)+
  theme_bw()+
  xlab("Myeloid response state") + ylab("") +
  theme(legend.position="none",
        axis.line = element_line(colour = "Black"),
        panel.background=element_rect(fill=NA),
        panel.border = element_rect(colour=NA),
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("Plots/Exhaustion.pdf", width=6.5,height=4.5,units="cm")



graphics.off()
rm(list = ls())


