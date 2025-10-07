BiocManager::install("GenomicFeatures")
BiocManager::install("ggplot2")
BiocManager::install("stringr")
BiocManager::install("tidyverse")
BiocManager::install("edgeR")

library("GenomicFeatures")
library("ggplot2")
library("stringr")
library("tidyverse")
library("edgeR")

setwd("C:/Users/User/Desktop/Università_Trento/Bioinformatics/Romanel/dataset")

## Data generated from GSE106305
## Load raw counts and gene annotation, stored in the "RNAseq_data.RData" file.
load("RNAseq_Data_PC3_OC2_Hypoxia.RData") 

## Check number of transcripts
dim(raw_counts_df)
head(raw_counts_df)
dim(r_anno_df)
head(r_anno_df)

## Check the library size of each sample
size_df <- data.frame("sample"=colnames(raw_counts_df),
                    "read_millions"=colSums(raw_counts_df)/1000000)

ggplot(data=size_df,aes(sample,read_millions)) +
   geom_bar(stat="identity",fill="grey50",colour="grey50", width=0.7, alpha=0.7) +
   coord_flip() +
   theme_bw()

## Create a data.frame with sample annotation (c_anno_df, separate conditions from replicates)
c_anno_df <- data.frame("sample"=colnames(raw_counts_df))
rownames(c_anno_df) <- c_anno_df$sample
c_anno_add <- as.data.frame(str_split(string=c_anno_df$sample, pattern="_", simplify = T))
colnames(c_anno_add) <- c("condition","replicate")
c_anno_df <- cbind(c_anno_df,c_anno_add)
c_anno_df

## Remove genes with low signal
## Filtered gene need to have at least #reads in at least #replicates in at least 1 condition

# count threshold
count_thr <- 20
# number of replicates with more counts than the count threshold
repl_thr <- 2 

filter_vec <- apply(raw_counts_df,1,
    function(y) max(by(y, c_anno_df$condition, function(x) sum(x>=count_thr))))
# see statistics for the filtering
table(filter_vec) 

filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,]
# check the dimension of the filtered matrix 
dim(filter_counts_df) 

# apply the filter on gene annotation
filter_anno_df <- r_anno_df[rownames(filter_counts_df),]
dim(filter_anno_df)

## Boxplot of gene counts after filtering absent genes

# gather is a function from the tidyr package
long_counts_df <- gather(filter_counts_df, key = "sample", value = "read_number")

# plot in log10 scale +1 to avoid infinite values
ggplot(data=long_counts_df,aes(sample,read_number+1)) +
   geom_boxplot(colour="indianred",fill="indianred",alpha=0.7) +
   theme_bw() +
   scale_y_log10()


## Check the library size of each sample after filtering
size_df <- data.frame("sample"=colnames(filter_counts_df),
                    "read_millions"=colSums(filter_counts_df)/1000000)

ggplot(data=size_df,aes(sample,read_millions)) +
   geom_bar(stat="identity",fill="indianred",colour="indianred",width=0.7,alpha=0.7)+
   coord_flip()+
   theme_bw()


## Cluster samples with hierarchical clustering and display the clustering tree

# scale for rows (each gene), to normalize for basal expression level differences
clu_data <- scale(t(filter_counts_df))
# distance matrix (we are clustering samples, not genes)
dd <- dist(clu_data, method = "euclidean") 
hc <- hclust(dd, method="ward.D")
# display dendogram
plot(hc) 

## PCA analysis of samples
data.matrix <- filter_counts_df
color <- rep("red",4)
color[which(c_anno_df$condition=="PC3siCtrlH")] = "blue"
data.PC <- prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1:2],xlim=c(-200,200),ylim=c(-200,200),col=color,pch=19)
text(data.PC$x[,1],data.PC$x[,2]+10,colnames(filter_counts_df),cex=0.7)

## DEG analysis with edgeR

# create a DGRList object
edge_c <- DGEList(counts=filter_counts_df,group=c_anno_df$condition,samples=c_anno_df,genes=filter_anno_df) 
edge_c

# normalization with the edgeR package (TMM method)
edge_n <- calcNormFactors(edge_c,method="TMM")
edge_n

# create a cpm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2)) #accounts for differences across samples, using a scaling factor
head(cpm_table)     #we take the matrix, calculate for each sample the cpm, at both inter and intra-samples level

# look at the boxplot distribution of gene expression signals after normalization
long_cpm_df <- gather(cpm_table, key = "sample", value = "CPM") 

ggplot(data=long_cpm_df,aes(sample,CPM+1)) +
   geom_boxplot(colour="olivedrab",fill="olivedrab",alpha=0.7)+
   theme_bw()+
   scale_y_log10() 


clu_data <- t(scale(t(cpm_table))) 
dd <- dist(t(clu_data), method = "euclidean") 
hc <- hclust(dd, method="ward.D")
plot(hc)    #difference across the two silenced samples

data.matrix <- cpm_table
color <- c(rep("red",2),rep("blue",2))
data.PC <- prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1:2],xlim=c(-200,200),ylim=c(-200,200),col=color,pch=19)
text(data.PC$x[,1],data.PC$x[,2]+10,colnames(filter_counts_df),cex=0.7)

# check ONECUT2 knockdown       #check what if the silencing was successful, extracting from the cpm table and create a barplot that cecks the expression of the silenced genes
onecut2 = cpm_table[which(rownames(cpm_table)=="ENSG00000119547"),]
long_onecut2_df <- gather(onecut2, key = "sample", value = "CPM") 

ggplot(data=long_onecut2_df,aes(sample,CPM+1)) +
   geom_bar(stat="identity",fill="indianred",colour="indianred",width=0.7,alpha=0.7)+
   theme_bw()

# define the experimental design matrix
design <- model.matrix(~0+group, data=edge_n$samples)     #same framework used for the array
colnames(design) <- levels(edge_n$samples$group)      #do a design model, creating a model.matrix that represents our design
rownames(design) <- edge_n$samples$sample     
design

# calculate dispersion and fit with edgeR (necessary for differential expression analysis)
edge_d <- estimateDisp(edge_n,design)   
edge_f <- glmQLFit(edge_d,design)   #perform the fitting

# definition of the contrast (conditions to be compared) 
contro <- makeContrasts("PC3siOC2H-PC3siCtrlH", levels=design)    #study variability and dispersion of the data

# fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)
DEGs <- as.data.frame(topTags(edge_t,n=20,p.value = 0.25,sort.by = "logFC"))  
DEGs <- as.data.frame(topTags(edge_t,n=20000))
DEGs$class <- "="
DEGs$class[which(DEGs$logCPM>1&DEGs$logFC>1.5&DEGs$FDR<0.25)] = "+"
DEGs$class[which(DEGs$logCPM>1&DEGs$logFC<(-1.5)&DEGs$FDR<0.25)] = "-"
DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]

head(DEGs)
table(DEGs$class)

## Display an MA plot
input_df <- DEGs
xlabel <- "log2 avg CPM (A)"
ylabel <- "log2 FC siOC2 vs siCtrl (M)"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$logCPM, input_df$logFC, xlab=xlabel, ylab=ylabel, 
     col=ifelse(input_df$class=="=","grey70","olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="MA plot")
abline(h=0,lty=2,col="grey20")

## Display a Volcano plot of the results:
input_df <- DEGs
xlabel <- "log2 FC siOC2 vs siCtrl"
ylabel <- "-log10 p-value"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$logFC, -log(input_df$PValue,base=10),xlab=xlabel, ylab=ylabel, 
     col=ifelse(input_df$class=="=","grey70","olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot")
abline(v=0,lty=2,col="grey20")

## Heatmap with DEG genes
cols <- c(rep("chartreuse4",2),rep("burlywood3",2)) 
pal <- c("blue","white","red") 
pal <- colorRampPalette(pal)(100)
heatmap(as.matrix(cpm_table[which(rownames(cpm_table)%in%DEGs$ensembl_gene_id[which(DEGs$class!="=")]),]),
        ColSideColors = cols,cexCol = 0.5,margins = c(4,4),col=pal,cexRow = 0.2)

## Export differentially expressed genes in a text file
up_DEGs <- DEGs[which(DEGs$class=="+"),]
down_DEGs <- DEGs[which(DEGs$class=="-"),]

write.table(up_DEGs,file="up_DEGs.txt",row.names=F,col.names=T,sep="/t",quote=F)
write.table(down_DEGs,file="down_DEGs.txt",row.names=F,col.names=T,sep="/t",quote=F)
write.table(DEGs,file="DEGs.txt",row.names=F,col.names=T,sep="/t",quote=F)
