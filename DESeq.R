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
library(openxlsx)
library(doBy) # for finding unique rows
library("ComplexHeatmap") # for heatmap


### Define directory for analysis ###
base.dir<-"/Volumes/qac/prj/csiszara/20181113-nmn/analysis/"

### Annotation files ###
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl") # link to biomart library for GRCh38 build
listAttributes(ensembl) # list of options to obtain from library
genes=getBM(attributes = c('ensembl_transcript_id','external_gene_name'),mart=ensembl) # selected options
colnames(genes)=c('Transcript ID', 'Gene') # name columns

########################
### Loading our data ###
########################
# point to kallisto files
kallisto.path<-paste0(base.dir,"03_kallisto") # path for all kallisto files
sample_id<-list.files(kallisto.path,pattern=".*",include.dirs=T,full.names=T) # list of sample IDs from the kallisto folders
s2c<-data.frame(path=as.character(sample_id),sample=basename(sample_id)) # data frame of kallisto sample info
s2c$path<-as.character(s2c$path) # structuring variable
my_read_tsv <- function(x, ...) { # function for reading tsv files
  read_tsv(x, guess_max = 40000, ...)
}
# Obtaining counts
txi <- tximport(paste0(s2c$path,"/abundance.tsv"), type = "kallisto", tx2gene = genes, countsFromAbundance="lengthScaledTPM")
raw.counts<-as.data.frame(txi$counts) # isolating raw counts
colnames(raw.counts)<-as.character(s2c$sample) # entering column names
# Keep raw counts
write.table(raw.counts,paste0(base.dir,"RawCounts.txt"), row.names=T, sep="\t", quote=F)


##############################################
### Creating subset with cuts for analysis ###
##############################################
# Remove counts that have all 0s for all samples or an overall average of 5 or less. This is removing "bad reads".
raw.counts.cut = raw.counts[-(which(rowSums(raw.counts) == 0 | rowMeans(raw.counts) <= 5)),]
raw.counts.cut = round(raw.counts.cut) # Rounding to integers for use in DESEq

#####################################################
### Comparison Day 30 vs Day 23 (case vs control) ###
#####################################################
condition<-data.frame(c(rep("AgedNMN",5), rep("Aged", 5),rep("Young",5)),  stringsAsFactors = default.stringsAsFactors()) # adding condition or treatment labels to samples
colnames(condition)<- 'Condition' # correctly labeling column for DESeq

### DESeq ###
dds<-DESeqDataSetFromMatrix(raw.counts.cut, # format data for DESeq analysis
                            condition,design=~ Condition) # matrix with counts, conditions or variables for the columns
# and design of experiment "~ 'variables included'" from column names of condition dataframe. Last variable listed is the ultimate comparison made.
dds <- DESeq(dds, betaPrior=FALSE) # Analysis, using DESeq dataset and not using zero-mean normal prior
dds_counts<-counts(dds,normalized=TRUE)  # collecting counts that are normalized by size factors
# write table of normalized counts
write.table(dds_counts,paste0(base.dir,"dds_counts.txt"),quote=F,row.names=T) 

#######################################
### PCA for further quality control ###
#######################################
groups = data.frame(c(rep("AgedNMN",5), rep("Aged", 5),rep("Young",5)),  stringsAsFactors = default.stringsAsFactors())
colnames(groups)<-"group"               #label column
####Start PCA plot ####
exp_data<-dds_counts #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_groups.png")) #set png destination
fviz_pca_ind(pca_exp_data, col.ind = groups$group , mean.point = FALSE, #Plot PCA with colored groups according to "groups" dataframe
             repel = TRUE)                                                                      
dev.off()  #png device off

#################################
### Obtain results from DESEq ###
#################################
### Comparisons data frames ###
comps = matrix(c('Aged','Young','AgedNMN','Young','AgedNMN','Aged'), nrow = 3, ncol = 2, byrow = TRUE) # create matrix of comparisons
comps = as.data.frame(comps)   # make a data frame
colnames(comps) = c('Group1', 'Group2')  # define column names

### Results ###
for(i in 1:nrow(comps)){     # For each comparison as in comps dataframe 
  res <- results(dds,cooksCutoff=F,contrast=c("Condition", paste0("",comps[i,1],""), paste0("",comps[i,2],""))) # obtain the results, comparison case then control
  res <-data.frame(res)      # Make results a dataframe
  Gene<-row.names(res)        # format row names
  res<-cbind(Gene,res)        # Add row names as column
  res$FC = 2^res$log2FoldChange # add fold change as a column
  res = cbind(res,dds_counts) # Add normalized counts to results file
  res = res[,c(1,9:23,2:3,8,4:7)] # rearrange columns so file in the format "ENSG, Gene, Normalized counts, basemean, log2FC, FC, lfcSE, stat,pvalue, adj p-value"
  res$DE<-ifelse(abs(res$log2FoldChange)>=log2(2) & res$padj<0.05,"TRUE","FALSE" )    # Detemine if meets DE criteria
  # write overall results
  write.xlsx(res,paste0(base.dir, "Results/DESeq2_",comps[i,1] ,"_",comps[i,2] ,".xlsx"))
  # write DE results
  write.xlsx(res[which(res$DE=="TRUE"),],paste0(base.dir,"Results/DESeq2_DE_",comps[i,1] ,"_",comps[i,2] ,".xlsx"))
} #end i


#### Further check on sample clustering ###

groups = data.frame(c(rep("AgedNMN",5), rep("Aged", 5),rep("Young",5)),  stringsAsFactors = default.stringsAsFactors())
colnames(groups)<-"group"               #label column
####Start PCA plot ####
exp_data<-dds_counts #read in counts
pca_exp_data<-prcomp(t(log(exp_data+1,2)),scale=F)  #run PCA
png(filename =paste0(base.dir,"Results/PCA_groups_2&3.png")) #set png destination
plot( pca_exp_data$x[,2], pca_exp_data$x[,3] )
text( pca_exp_data$x[,2], pca_exp_data$x[,3], colnames( exp_data ) )                                                               
dev.off()  #png device off

png(filename =paste0(base.dir,"Results/Screeplot.png")) #set png destination
fviz_eig(pca_exp_data)
dev.off()

############################################################################################
### Heatmap of 100 from each analysis with pseudogenes removed and displaying raw counts ###
############################################################################################
genes2=getBM(attributes = c('ensembl_transcript_id','external_gene_name', 'description'),mart=ensembl) # selected options
colnames(genes2)=c('Transcript ID', 'Gene',"Description") # name columns
### Read in results and format ###
AY=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_Aged_Young.xlsx"))
AY = AY[order(AY$padj),] # Order results by padj 
AY$Description = genes2[match(AY$Gene , genes2$Gene),3] # annotate results with description
AY = AY[-grep("pseudogene", AY$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
AY = AY[1:100,] # Obtain top 100 DE genes from these results
NMNA=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_AgedNMN_Aged.xlsx"))
NMNA = NMNA[order(NMNA$padj),] # Order results by padj 
NMNA$Description = genes2[match(NMNA$Gene , genes2$Gene),3] # annotate results with description
NMNA = NMNA[-grep("pseudogene", NMNA$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
NMNA = NMNA[1:100,] # Obtain top 100 DE genes from these results
NMNY=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_AgedNMN_Young.xlsx"))
NMNY = NMNY[order(NMNY$padj),] # Order results by padj 
NMNY$Description = genes2[match(NMNY$Gene , genes2$Gene),3] # annotate results with description
NMNY = NMNY[-grep("pseudogene", NMNY$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
### merging sets together
top = rbind(AY,NMNA,NMNY) # combine datasets
top = top[firstobs(top$Gene),] # remove duplicates
rownames(top)=top$Gene
top = top[,-1]
colnames(top)=c("AgedNMN1.0", "AgedNMN1.1", "AgedNMN2.0", 'AgedNMN3.1', "AgedNMN5.0", "Aged1.2", "Aged2.1", "Aged3.1","Aged8.1",
                "Aged9.1", "Young1.0", "Young1.1", "Young1.2", "Young1.4", "Young2.0")
top = as.matrix(top)
### Heatmap
png(paste0(base.dir,"Results/Top100_Heatmap.png"),height=2024, width=1764)
heatmap(top,cexCol=2.5, cexRow=1)
dev.off()

#####################################################################################################
### Heatmap of all DE genes from each analysis with pseudogenes removed and displaying raw counts ###
#####################################################################################################
### Read in results and format ###
AY=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_Aged_Young.xlsx"))
AY$Description = genes2[match(AY$Gene , genes2$Gene),3] # annotate results with description
AY = AY[-grep("pseudogene", AY$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
NMNA=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_AgedNMN_Aged.xlsx"))
NMNA$Description = genes2[match(NMNA$Gene , genes2$Gene),3] # annotate results with description
NMNA = NMNA[-grep("pseudogene", NMNA$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
NMNY=read.xlsx(paste0(base.dir,"Results/DESeq2_DE_AgedNMN_Young.xlsx"))
NMNY$Description = genes2[match(NMNY$Gene , genes2$Gene),3] # annotate results with description
NMNY = NMNY[-grep("pseudogene", NMNY$Description),c(1:16)]  # remove genes that are pseudogenes and prune to only genes and counts
### merging sets together
top = rbind(AY,NMNA,NMNY) # combine datasets
top = top[firstobs(top$Gene),] # remove duplicates
rownames(top)=top$Gene
top = top[,-1]
colnames(top)=c("AgedNMN1.0", "AgedNMN1.1", "AgedNMN2.0", 'AgedNMN3.1', "AgedNMN5.0", "Aged1.2", "Aged2.1", "Aged3.1","Aged8.1",
                "Aged9.1", "Young1.0", "Young1.1", "Young1.2", "Young1.4", "Young2.0")
top = as.matrix(top)
### Heatmap
png(paste0(base.dir,"Results/Heatmap.png"),height=2024, width=1764)
heatmap(top,cexCol=2.5, cexRow=1)
dev.off()



