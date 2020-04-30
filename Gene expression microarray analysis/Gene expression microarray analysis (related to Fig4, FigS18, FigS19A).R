#############################
######### Load data #########
#############################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("limma","biomaRt","ExpressionNormalizationWorkflow","sva","snm", "Biobase","AnnotationDbi","ggplot2","clusterProfiler","GSEABase","enrichplot","org.Hs.eg.db","GOSemSim","RColorBrewer","cowplot","magrittr"))

library(limma)
## Set your working directory (or use Control + Shift + H to choose it) ##
## The working directory should contain the RData file "Gene expression microarray analysis.RData" ###
setwd("...")

### load data
load("Gene expression microarray analysis.RData")

## check if you properly have 62976 probes and 21 samples ##
dim(x)

#############################
###### Gene Anotation #######
#############################

## check if you properly have 62976 probe annotation ##
length(unique(annot.info$ID))

## merge expression data and annotation info ##
x$genes <- merge(x$genes, annot.info, by.x = "FeatureNum", 
                 by.y = "ID", sort = FALSE)


#############################
### Background correction ### 
#############################

y <- backgroundCorrect(x, method="normexp")
y$E <- log2(y$E+1)
boxplot(y$E,main="expression value",las=2)
#View(y$genes)


#############################
####### Gene filtering ###### 
#############################

## Gene filtering ##
isControl <- y$genes$CONTROL_TYPE!= "FALSE"
yfilt <- y[!isControl, ]
dim(yfilt)

## clean up the SampleNames ##
colnames(yfilt$E) <- sub("-T@.*", "", colnames(yfilt$E))

## combine data measured by probes with same ProbeIDs ##
UniProbe <- aggregate(x = yfilt$E, 
                      by = list(yfilt$genes$ProbeName), FUN = mean)
colnames(UniProbe)[1] <- "ProbeNames"
#View(UniProbe)
dim(UniProbe)
boxplot(UniProbe[,-1],las=2)

#############################
####### Normalization #######
#############################

### First, evaluate the batch effects ###
#########################################

library(ExpressionNormalizationWorkflow)
library(sva)
library(Biobase)

## Create an Expression Set object ##
# Re-arrang expression data
exprs <- UniProbe
row.names(exprs) <- exprs[,1]
exprs[,1] <- NULL
#View(exprsAll)

# import patient=(sample) information
#covrts <-  read.csv("HCC sample info.csv",
#                    header=TRUE, as.is=TRUE, row.names=1)

# re-order the expression matrix and the sample matrix
exprs <- exprs[,order(colnames(exprs))]
covrts <- covrts[order(row.names(covrts)),]

# construct ExpressionSet object
inpData <- expSetobj(exprs, covrts)
annotation(inpData) <- "org.Hs.eg.db"
dim(inpData)

## Set the covariates whose effect size on the data needs to be calculated (optional)
cvrts_eff_var <- c("Batch","Group","RnaQC","GoodPrognosis",
                   "No.Sex","PDL1.Mf","PDL1.TC")
## Set a PVCA Threshold Value between 0 & 1
## PVCA Threshold Value is the percentile value of the minimum amount of the variabilities that the selected principal components need to explain, here requiring 75% of the expression variance to be captured by the PCs
pct_thrsh <- 0.75
## Perform the PVCA analysis
pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)


### use "sva" method to removing batch   ###
### effects and other unwanted variation ###
############################################

### Setting up the data ###
pheno = pData(inpData)
edata = exprs(inpData)

mod = model.matrix(~as.factor(pheno$Group)+
                    as.factor(Batch) +
                    as.factor(No.Sex) +
                    as.factor(pheno$PDL1.Mf) +
                    as.factor(pheno$PDL1.TC),data=pheno)
mod0 = model.matrix(~as.factor(Batch) +
                     as.factor(No.Sex) +
                     as.factor(pheno$PDL1.Mf) +
                     as.factor(pheno$PDL1.TC),data=pheno)

### Applying the sva function to  ########
### estimate batch and other artifacts ###
svobj = sva(edata,mod,mod0)
colnames(svobj$sv) <- paste("sv",1:svobj$n.sv,sep="")
phenoSv <- cbind(pheno,svobj$sv)


### Computing the correlation between the  ###
### surrogate variables and the covariates ###
##############################################

glm.sv1 <- glm(phenoSv[,"sv1"]~
                 phenoSv[,"Group"] + phenoSv[,"GoodPrognosis"]+
                 phenoSv[,"Batch"] + phenoSv[,"RnaQC"] +
                 phenoSv[,"No.Sex"] + phenoSv[,"PDL1.Mf"] +
                 phenoSv[,"PDL1.TC"])
summary(glm.sv1)

glm.sv2 <- glm(phenoSv[,"sv2"]~
                 phenoSv[,"Group"] + phenoSv[,"GoodPrognosis"]+
                 phenoSv[,"Batch"] + phenoSv[,"RnaQC"] +
                 phenoSv[,"No.Sex"] + phenoSv[,"PDL1.Mf"] +
                 phenoSv[,"PDL1.TC"]) 
summary(glm.sv2)

glm.sv3 <- glm(phenoSv[,"sv3"]~
                 phenoSv[,"Group"] + phenoSv[,"GoodPrognosis"]+
                 phenoSv[,"Batch"] + phenoSv[,"RnaQC"] +
                 phenoSv[,"No.Sex"] + phenoSv[,"PDL1.Mf"] +
                 phenoSv[,"PDL1.TC"]) 
summary(glm.sv3)

glm.sv4 <- glm(phenoSv[,"sv4"]~
                 phenoSv[,"Group"] + phenoSv[,"GoodPrognosis"]+
                 phenoSv[,"Batch"] + phenoSv[,"RnaQC"] +
                 phenoSv[,"No.Sex"] + phenoSv[,"PDL1.Mf"] +
                 phenoSv[,"PDL1.TC"]) 
summary(glm.sv4)

glm.sv5 <- glm(phenoSv[,"sv5"]~
                 phenoSv[,"Group"] + phenoSv[,"GoodPrognosis"]+
                 phenoSv[,"Batch"] + phenoSv[,"RnaQC"] +
                 phenoSv[,"No.Sex"] + phenoSv[,"PDL1.Mf"] +
                 phenoSv[,"PDL1.TC"]) 
summary(glm.sv5)


### The results show that sv1 is correlated with Batch
### sv2 is correlated with PDL1.Mf
### sv3 is weakly correlated with sex. 
############################################################



### Supervised normalization of Microarrays ###
###############################################

### Use independent SVs as adjustment variable
## Choose the biological variableof interest
bv <- c("Group")
## Choose your adjustment variable of interest, 
#  combine 'Batch' and independent SVs
av <- c("sv1","sv3","sv4","sv5") 
## The intensity-dependent adjustment variables adjust for array effects 
iv <- c("Array") 
## Run SNM
sv_snmObj <- snmAnaly(edata, phenoSv, bv, av, iv)
## After comparison with other models, we chose this one for normalization
graphics.off()


## Create an expressionSet object of the normalized dataset(s)
sv_snmNorm_data <- sv_snmObj$norm.dat
#View(phenoSv)
dim(sv_snmNorm_data)
colnames(sv_snmNorm_data) <- rownames(phenoSv)
#View(sv_snmNorm_data)
sv_snm_data <- expSetobj(sv_snmNorm_data, phenoSv)
boxplot(sv_snmNorm_data,las=2)

## Write this dateset to a table with rows as genes and columns as samples (with names the same as that from the input file)
## By doing this, you may pause the analaysis here
write.csv(sv_snmNorm_data, file="NormalizedData.csv")
write.csv(phenoSv, file="Sample info for NormalizedData.csv")


#####################################################
### (optional) Annotate data for other analayses ####
#####################################################

#############################################################
### The following codes can help you Skip above precedure ###
## Set your working directory (or use Control + Shift + H to choose it)
#setwd("...")

## Read normalized expression data
sv_snmNorm_data <- read.csv("NormalizedData.csv",
                           header=TRUE, as.is=TRUE, row.names=1)

## select useful annotation info
UseAnnotInfo <- annot.info[,c("NAME","UNIGENE_ID")]
dim(UseAnnotInfo)
UseAnnotInfo <- UseAnnotInfo[!duplicated(UseAnnotInfo),]
dim(UseAnnotInfo)

## annotate normalized data (transfer probe names into unigene IDs)
sv_snmNorm_data <- as.data.frame(sv_snmNorm_data)
sv_snmNorm_data$NAME <- rownames(sv_snmNorm_data)
Norm_data.Annot <- merge(sv_snmNorm_data, UseAnnotInfo, 
                         all.y = FALSE,
                         all.x=FALSE,sort = FALSE)
sv_snmNorm_data$NAME <- NULL
Norm_data.Annot$NAME <- NULL
Norm_data.Annot <- Norm_data.Annot[!is.na(Norm_data.Annot$UNIGENE_ID),]
Norm_data.Annot <- aggregate(x = Norm_data.Annot[,-ncol(Norm_data.Annot)], 
                                by = list(Norm_data.Annot[,ncol(Norm_data.Annot)]), 
                                FUN = median)
colnames(Norm_data.Annot)[1] <- "UNIGENE_ID"
dim(Norm_data.Annot)
#View(Norm_data.Annot)

## (optional) annotate normalized data for CIBERSORT
library(AnnotationDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
CIBER_data.Annot <- Norm_data.Annot
CIBER_data.Annot[,1] <- mapIds(org.Hs.eg.db,
                              keys=as.character(CIBER_data.Annot$UNIGENE_ID),
                              column="SYMBOL",
                              keytype="UNIGENE",
                              multiVals="first")
colnames(CIBER_data.Annot)[1] <- "SYMBOL"
dim(CIBER_data.Annot)
CIBER_data.UniAnnot <- CIBER_data.Annot[!is.na(CIBER_data.Annot[,1]),]
dim(CIBER_data.UniAnnot)
CIBER_data.UniAnnot <- aggregate(x = CIBER_data.Annot[,-1], 
                        by = list(CIBER_data.Annot[,1]), FUN = median)
dim(CIBER_data.UniAnnot)
dim(CIBER_data.UniAnnot)
## Export the annotated data ##
write.table(CIBER_data.UniAnnot,
            file="DataForCIBERSORT.txt",
            sep="\t", col.names=T, row.names=F,quote=FALSE)



## (optional) export normalized data annotated using ExtrezGene ID
library(AnnotationDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
Entrez_data.Annot <- Norm_data.Annot
Entrez_data.Annot[,1] <- mapIds(org.Hs.eg.db,
                              keys=as.character(Entrez_data.Annot$UNIGENE_ID),
                              column="ENTREZID",
                              keytype="UNIGENE",
                              multiVals="first")
colnames(Entrez_data.Annot)[1] <- "ENTREZID"
dim(Entrez_data.Annot)
Entrez_data.UniAnnot <- Entrez_data.Annot[!is.na(Entrez_data.Annot[,1]),]
dim(Entrez_data.UniAnnot)
Entrez_data.UniAnnot <- aggregate(x = Entrez_data.Annot[,-1], 
                                  by = list(Entrez_data.Annot[,1]), 
                                  FUN = median)
dim(Entrez_data.UniAnnot)
rownames(Entrez_data.UniAnnot) <- Entrez_data.UniAnnot[,1]
dim(Entrez_data.UniAnnot)
## Export the annotated data ##
write.csv(Entrez_data.UniAnnot[,-1], row.names=TRUE,
          file="NormallizedData-EntrezGene.csv")


##########################################
### Differential Expresssion Anlaysis ####
##########################################

## Load necessary R packages for this section
library(ExpressionNormalizationWorkflow)
library(Biobase)
library(limma)

## Set your working directory (or use Control + Shift + H to choose it)
#setwd("...")

## Read phenotype data
Pheno <-  read.csv("Sample info for NormalizedData.csv",
                   header=TRUE, as.is=TRUE, row.names=1)

## Read normalized expression data
sv_snmNorm_data <- read.csv("NormalizedData.csv",
                           header=TRUE, as.is=TRUE, row.names=1)



### Construct ExpresssionSet ###
sv_snm_data <- expSetobj(sv_snmNorm_data, Pheno)
annotation(sv_snm_data) <- "org.Hs.eg.db"

### Differential Expresssion Anlaysis ###
design1 <- model.matrix(~0+as.factor(Group),data = pData(sv_snm_data))
colnames(design1) <- c("low","hi")
fit1 <- lmFit(exprs(sv_snm_data), design1)

cont.matrix1 <- makeContrasts("HIvsLOW"=hi-low, levels=design1)
cfit1 <- contrasts.fit(fit1, cont.matrix1)
efit1 <- eBayes(cfit1, trend=TRUE, robust=TRUE)
results1 <- decideTests(efit1)
summary(results1)
Output1 <- topTable(efit1, coef=1, number=Inf)
dim(Output1)


#########################
### Annotate Results ####
#########################

Output1$NAME <- row.names(Output1)
Output1$Reg <- Output1$logFC >0
Output1$Reg <- sub("TRUE","Up", Output1$Reg)
Output1$Reg <- sub("FALSE","Down", Output1$Reg)
Output1$FC <- 2^Output1$logFC
#View(Output1)


## Read probes annotation info
UseAnnotInfo <- annot.info[,c("NAME","UNIGENE_ID")]
dim(UseAnnotInfo)
UseAnnotInfo <- UseAnnotInfo[!duplicated(UseAnnotInfo),]
dim(UseAnnotInfo)

## annotate differentially expressing probes 
## (translate probe names into unigene IDs)
PrAnnotResult1 <- merge(Output1, UseAnnotInfo, all.y = FALSE,
                       sort = FALSE)
dim(PrAnnotResult1)
colnames(PrAnnotResult1)[1] <- "ProbeID"

## exclude differentially expressing probes without Unigene IDs
PrAnnotResult1 <- PrAnnotResult1[!is.na(PrAnnotResult1$UNIGENE_ID),]
length(PrAnnotResult1$UNIGENE_ID)
length(unique(PrAnnotResult1$UNIGENE_ID))

## annotate genes with EntrezGene IDs
library(AnnotationDbi)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
EnAnnotResult1 <- PrAnnotResult1

EnAnnotResult1$EntrezID <- mapIds(org.Hs.eg.db,
                                 keys=as.character(EnAnnotResult1$UNIGENE_ID),
                                 column="ENTREZID",
                                 keytype="UNIGENE",
                                 multiVals="first")
EnAnnotResult1$Symbol <- mapIds(org.Hs.eg.db,
                               keys=as.character(EnAnnotResult1$UNIGENE_ID),
                               column="SYMBOL",
                               keytype="UNIGENE",
                               multiVals="first")
EnAnnotResult1$GeneName <- mapIds(org.Hs.eg.db,
                                 keys=as.character(EnAnnotResult1$UNIGENE_ID),
                                 column="GENENAME",
                                 keytype="UNIGENE",
                                 multiVals="first")

## exclude differentially expressing probes without EntrezGene IDs
EnAnnotResult1 <- EnAnnotResult1[!is.na(EnAnnotResult1$EntrezID),]
#View(EnAnnotResult1)
dim(EnAnnotResult1)


## keep only one result data related to each EntrezID
EnAnnotResult1 <- EnAnnotResult1[order(EnAnnotResult1$adj.P.Val),]
EnAnnotResult1 <- EnAnnotResult1[!duplicated(EnAnnotResult1$EntrezID),]
dim(EnAnnotResult1)

## re-order columns
OutputGenes1 <- cbind(EnAnnotResult1$EntrezID,EnAnnotResult1$UNIGENE_ID,
                      EnAnnotResult1$Symbol,EnAnnotResult1$GeneName,
                      EnAnnotResult1$ProbeID,EnAnnotResult1[2:8])
colnames(OutputGenes1) <- gsub("EnAnnotResult1\\$", "", colnames(OutputGenes1))
#View(OutputGenes1)

### Export the result of the differential expression analysis ###
write.csv(OutputGenes1, file="Differential expression analysis result.csv",
          row.names = FALSE)

### Pause point. You may run remove(list = ls()) to clean the runing environment. ###
#remove(list = ls())

##################
### load DEGs ####
##################

### load Expression Microarray result data 
#setwd("...")
Genes <- read.csv("Differential expression analysis result.csv",
                   header=TRUE, as.is=TRUE)

######################
### Volcanol plot ####
######################
library(ggplot2)
Genes$Group <- NA
Genes$Group <- ifelse(Genes$adj.P.Val<0.05 & abs(Genes$logFC)>1,"Sig","notSig")
Genes$Group[Genes$adj.P.Val<0.05 & abs(Genes$logFC)>1] <- ifelse(Genes$logFC[Genes$adj.P.Val<0.05 & abs(Genes$logFC)>1]>0,"Up","Dn")
Genes$Group[Genes$adj.P.Val<0.01 & abs(Genes$logFC)>1] <- ifelse(Genes$logFC[Genes$adj.P.Val<0.01 & abs(Genes$logFC)>1]>0,"Upstr","Dnstr")
Genes$Group <- factor(Genes$Group, levels=c("notSig","Dn","Up","Dnstr","Upstr"))
table(Genes$Group)

ggplot(data=Genes, aes(x=logFC, y=-log10(adj.P.Val), colour=Group,size=Group)) +
  geom_point(alpha=.4) +
  scale_size_manual(values=c(.2,.2,.5,.2,.5)) +
  scale_color_manual(values=c("grey70","darkblue","darkred","#0571b0","#ca0020")) +
  xlab("log2 Fold Change") + ylab("-log10 P-value")+
  theme(legend.position="none",
        axis.title = element_text(size=8,colour="black"),
        axis.text = element_text(size=8,colour="black"))
ggsave("DEG volcano.pdf", width=7,height=7,units="cm")
dev.off()

###############################
### Gene ontology analyses ####
###############################

### select genes with ajusted p values <= 0.05 & logFC>=1
degGenes <- Genes[Genes$adj.P.Val<=0.05 & abs(Genes$logFC)>=1,]


### construct the gene list ####
### select genes with ajusted p values <= 0.05 & logFC>=1
GenesGL <- 2^(degGenes$logFC)
names(GenesGL) <- as.character(degGenes$EntrezID)
GenesGL <- sort(GenesGL, decreasing = TRUE)
length(GenesGL)


## load GO-biological processes gene sets
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(magrittr)
require(RColorBrewer)

## GO analysis - Up regulated
egoGenesUP <- enrichGO(names(GenesGL[GenesGL>1]), OrgDb = "org.Hs.eg.db", 
                        ont="BP", readable=TRUE)
goplot(egoGenesUP)
egoGenesUP@result$SetSize <- as.numeric(gsub("/.*", "",egoGenesUP@result$BgRatio))
egoGenesUP@result$FoldEnrich <- egoGenesUP@result$Count/egoGenesUP@result$SetSize
#View(egoGenesUP@result)
write.table(egoGenesUP@result,
            file="GO analysis results.txt",
            sep="\t", col.names=T, row.names=F,quote=FALSE)



## Dotplot of the Up-GO analysis result
MyPal <- c(brewer.pal(12,"Set3")[c(3,6,8,10,11)],"#A6CEE3",brewer.pal(8,"Set2")[-6])
brewer.pal(9,"YlOrRd")
egoUpPlot <- egoGenesUP@result[1:25,]
MRSupDot <- ggplot(egoUpPlot,aes(x=-log(p.adjust),
                                 y=reorder(Description, -log(p.adjust)),
                                 size=Count,colour= (FoldEnrich))) +
                  geom_point()+
                  scale_size_continuous(range=c(2,6))+
                  scale_color_gradientn(colours=brewer.pal(9,"YlOrRd")[3:9])+
                  theme_bw()
MRSupDot

## cnetplot of the Up-GO analysis result
#graphics.off()
cnetplot(egoGenesUP,showCategory = 25,fixed=FALSE,
         foldChange=NULL)
egoGenesUPImm <- dropGO(egoGenesUP,term=c("GO:0030198"))
cnetplot(egoGenesUPImm,showCategory = 24,fixed=FALSE,foldChange=NULL) +
  theme(legend.position="none") +
  coord_flip() +
  scale_y_reverse()

## GO analysis - Down regulated
egoGenesDN <- enrichGO(names(GenesGL[GenesGL>1]), OrgDb = "org.Hs.eg.db", 
                        ont="BP", readable=TRUE)
goplot(egoGenesDN)
egoGenesDN@result$SetSize <- as.numeric(gsub("/.*", "",egoGenesDN@result$BgRatio))
egoGenesDN@result$FoldEnrich <- egoGenesDN@result$Count/egoGenesDN@result$SetSize
View(egoGenesDN@result)

## Dotplot of the DN-GO analysis result
egoDNPlot <- egoGenesDN@result[1:25,]
MRSDNDot <- ggplot(egoDNPlot,aes(x=-log(p.adjust),
                                 y=reorder(Description, -log(p.adjust)),
                                 size=Count,colour= (FoldEnrich))) +
                    geom_point()+
                    scale_size_continuous(range=c(2,6))+
                    scale_color_gradientn(colours=brewer.pal(9,"YlOrRd")[3:9])+
                    theme_bw()
MRSDNDot


#####################
### GSE analyses ####
#####################

grpGL <- Genes$logFC 
names(grpGL) <- as.character(Genes$EntrezID)
grpGL <- sort(grpGL, decreasing = TRUE)
length(grpGL)

## Immune gene sets GESA
# ImmGS <-  read.delim("Gene sets collection/Immune gene sets.txt",
#                    header=TRUE, as.is=TRUE,sep = "\t")
ImmGS <- ImmGS[,c(1,3)]
ImmGSEA <- GSEA(grpGL, TERM2GENE=ImmGS, minGSSize=1,
               maxGSSize=Inf,verbose=FALSE)
ImmGSEA <- setReadable(ImmGSEA, OrgDb = org.Hs.eg.db,keytype="ENTREZID")
View(ImmGSEA@result)

#Exhausted T cells
#NES = 2.383568, p = 0.002542373
gseaplot(ImmGSEA, 1)

#TAM
#NES = 2.191028, p = 0.002542373
gseaplot(ImmGSEA, 2)

#MDSC
#NES = 1.816189, p = 0.006861063
gseaplot(ImmGSEA, 3)


graphics.off()
rm(list = ls())