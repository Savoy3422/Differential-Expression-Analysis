##############################
### Load required packages ###
##############################
require("sleuth")
require("tximport")
library("readr")
library("IRanges")
require("DESeq2")
library(biomaRt)
library(factoextra)
library(pca3d) # for 3d PCA
library(openxlsx)


### Define directory for analysis ###
base.dir<-"/Volumes/qac/prj/csiszara/20181113-hp/analysis/"

### Annotation files ###
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl") # link to biomart library for GRCh38 build
#listAttributes(ensembl) # list of options to obtain from library
genes=getBM(attributes = c('ensembl_gene_id','external_gene_name', 'description'),mart=ensembl) # selected options
colnames(genes)=c('Ensembl_ID', 'Gene', "Description") # name columns

annot<-read.table("/Volumes/acbm/data/Indra_documents/STAR_GRCh38_reference/Gene_Mus_musculus.GRCm38.88.clean.txt",header=T,sep=" ")


########################
### Loading our data ###
########################
### Bluebee Lexogen pipeline used to obtain counts ###
quantseq_samples<-read.table(paste0(base.dir,"Bluebee/sampleid.txt"),sep="\t",header=F,row.names = 1) # obtain sample names
raw.counts<-matrix(0,nrow=47648,ncol=nrow(quantseq_samples)) # nrow is the total number of rows in the read_couts.txt ( wc -l read_counts.txt) read counts file. ncol is the number of samples
for (i in 1:nrow(quantseq_samples)){ # for each sample
  sample_id<-row.names(quantseq_samples)[i] # pick out the sample id
  dir_id<-paste(base.dir,"Bluebee/",sample_id,"/star_out/", sample_id, "_R1_001.fastq.gz/read_counts.txt",sep="") # create path for counts
  rnaSeq<-read.table(dir_id,header=F) # read counts data
  raw.counts[,i]<-rnaSeq[,2] # putting counts into a column for the sample
}
raw.counts = as.data.frame(raw.counts) # change to dataframe
raw.counts = cbind(rnaSeq$V1,raw.counts) # add ENSG numbers
colnames(raw.counts) = c("Ensemble_gene", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                         "16", "17", "18", "19", "20")
row.names(raw.counts) = raw.counts$Ensemble_gene # changing row names for recall after DESeq
raw.counts$Gene = genes$Gene[match(rownames(raw.counts), genes$Ensembl_ID)] # annotate genes based on ENSG number
raw.counts = raw.counts[,c(1,ncol(raw.counts),2:(ncol(raw.counts)-1))] # rearrange columns
# save raw counts
write.table(raw.counts, paste0(base.dir,"rawCounts.txt"),row.names=F, sep="\t", quote=F)


#######################################################
### Creating groups labels to be compared for DESeq ###
#######################################################
counts=raw.counts[-(which(rowMeans(raw.counts[,3:ncol(raw.counts)])<= 10 | rowSums(raw.counts[,3:ncol(raw.counts)])== 0)),] # remove low reads from analysis
condition<-data.frame(c(rep("Young",5), rep("Aged",5), rep("YoungAged",5), rep("AgedYoung",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
counts<-counts[,3:ncol(counts)] # trimming raw counts for DESeq (genes not needed)

#############
### DESeq ###
#############
dds<-DESeqDataSetFromMatrix(counts, # format data for DESeq analysis
                            condition,design=~condition) # matrix with counts, conditions or variables for the columns
# and design of experiment "~ 'variables included'"
dds <- DESeq(dds, betaPrior=FALSE) # Analysis, using DESeq dataset and not using zero-mean normal prior
dds_counts<-counts(dds,normalized=TRUE)  # collecting counts that are normalized by size factors
# write table of normalized counts
write.table(dds_counts,paste0(base.dir,"dds_counts.txt"),quote=F,row.names=T) 

#######################################
### PCA for further quality control ###
#######################################
####Start PCA plot of all groups ####
exp_data<-dds_counts #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_all_groups.png"), width=764, height=764) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

####Start PCA plot of all groups 2 and 9 removed ####
condition<-data.frame(c(rep("Young",4), rep("Aged",4), rep("YoungAged",5), rep("AgedYoung",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,-c(2,9)] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_all_groups_2_9_removed.png"), width=764, height=764) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start 3D PCA Plot with all groups ###
cols = c(rep("blue",5), rep("red",5), rep("purple",5), rep("green",5))
exp_data = dds_counts
pca_exp_data<-prcomp(t(log2(exp_data+1)),scale=F)
pca3d(pca_exp_data, show.labels = TRUE)
snapshotPCA3d(file=paste0(base.dir,"Results/PCA_3d.png"))

#### Start PCA plot of only Young and Aged ####
condition<-data.frame(c(rep("Young",5), rep("Aged",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,1:10] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_Young_and_Aged_groups.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start PCA plot of only YoungAged and AgedYoung ####
condition<-data.frame(c(rep("YoungAged",5), rep("AgedYoung",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,11:20] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_YoungAged_and_AgedYoung_groups.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start PCA plot of only Young and Aged removing samples 1 and 7 ####
condition<-data.frame(c(rep("Young",4), rep("Aged",4))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,c(2:6, 8:10)] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_Young_and_Aged_groups_1_7_removed.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start PCA plot of only Young and Aged removing samples 2 and 9 ####
condition<-data.frame(c(rep("Young",4), rep("Aged",4))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,c(1,3:8,10)] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_Young_and_Aged_groups_2_9_removed.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start PCA plot of only Young, Aged, and YoungAged ####
condition<-data.frame(c(rep("Young",5), rep("Aged",5), rep("YoungAged",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,1:15] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_Young_Aged_and_YoungAged_groups.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#### Start PCA plot of only Young, Aged, and AgedYoung####
condition<-data.frame(c(rep("Young",5), rep("Aged",5), rep("AgedYoung",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable
exp_data<-dds_counts[,c(1:10,16:20)] #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_Young_Aged_and_AgedYoung_groups.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = condition$condition , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off


#######################################################
### Creating groups labels to be compared for DESeq ###
#######################################################
counts<-counts[,-(c(2,9))] # pruning out samples 2 and 9 
condition<-data.frame(c(rep("Young",4), rep("Aged",4), rep("YoungAged",5), rep("AgedYoung",5))) # adding condition or treatment labels to samples
colnames(condition)= "condition" # structuring variable


#############
### DESeq ###
#############
dds<-DESeqDataSetFromMatrix(counts, # format data for DESeq analysis
                            condition,design=~condition) # matrix with counts, conditions or variables for the columns
# and design of experiment "~ 'variables included'"
# filtering genes not expressed uniformly in either group, or study-wide 
goodRows<-which(apply(counts(dds),1,function(x) (all(x[which(condition$condition=="Young")]>0) & all(x[which(condition$condition=="Aged")]>0)) & 
                        all(x[which(condition$condition=="AgedYoung")]>0) & all(x[which(condition$condition=="YoungAged")]>0))  & rowSums(counts(dds)) >= 10)
dds<-dds[goodRows,]
dds <- DESeq(dds, betaPrior=FALSE) # Analysis, using DESeq dataset and not using zero-mean normal prior
dds_counts<-counts(dds,normalized=TRUE)  # collecting counts that are normalized by size factors
# write table of normalized counts
write.table(dds_counts,paste0(base.dir,"dds_counts_no2_9.txt"),quote=F,row.names=T) 

### Comparisons data frames ###
comps = matrix(c('YoungAged','Young','AgedYoung','Aged','Aged','Young'), nrow = 3, ncol = 2, byrow = TRUE) # create matrix of comparisons
comps = as.data.frame(comps)   # make a data frame
colnames(comps) = c('Group1', 'Group2')  # define column names

### Obtain comparison ALL results from DESeq analysis ###
for(i in 1:nrow(comps)){     # For each comparison as in comps dataframe 
  res <- results(dds,cooksCutoff=F,contrast=c("condition", paste0("",comps[i,1],""), paste0("",comps[i,2],""))) # obtain the results, comparison case then control
  res <-data.frame(res)      # Make results a dataframe
  Ensembl_ID<-row.names(res)        # format row names
  res<-cbind(Ensembl_ID,res)        # Add row names as column
  res$FC = 2^res$log2FoldChange # add fold change as a column
  res = cbind(res,dds_counts) # Add normalized counts to results file
  res$Gene = genes[match(rownames(res) , genes$Ensembl_ID),2] # annotate results with genes
  res = res[,c(1,27,9:26,2:3,8,4:7)] # rearrange columns so file in the format "ENSG, Gene, Normalized counts, basemean, log2FC, FC, lfcSE, stat,pvalue, adj p-value"
  res$DE<-ifelse(abs(res$log2FoldChange)>=log2(2) & res$padj<0.05,"TRUE","FALSE" )    # Detemine if meets DE criteria
  # write overall results
  write.xlsx(res,paste0(base.dir, "Results/DESeq2_",comps[i,1] ,"_vs_",comps[i,2] ,".xlsx"))
  # write DE results
  write.xlsx(res[which(res$DE=="TRUE"),],paste0(base.dir,"Results/DESeq2_DE_",comps[i,1] ,"_vs_",comps[i,2] ,".xlsx"))
} #end i


##################
### Analysis A ###
##################
 # No overlap of genes

##################
### Analysis B ###
##################
DE3 = read.xlsx(paste0(base.dir,"Results/DESeq2_DE_Aged_vs_Young.xlsx"))
row.names(DE3) = DE3$Gene
DE3 = DE3[,c(3:20)]
colnames(DE3) = c("Young1", "Young3", "Young4", "Young5", "Aged6", "Aged7", "Aged8", "Aged10",
                  "YoungAged11", "YoungAged12", "YoungAged13", "YoungAged14", "YoungAged15", "AgedYoung16",
                  "AgedYoung17", "AgedYoung18", "AgedYoung19", "AgedYoung20")
DE3 = as.matrix(DE3)
png(paste0(base.dir,"Results/Aged_vs_Young_DE_Heatmap_Counts.png"),height=2024, width=1764)
heatmap(DE3,cexCol=2.5, cexRow=1)
dev.off()

