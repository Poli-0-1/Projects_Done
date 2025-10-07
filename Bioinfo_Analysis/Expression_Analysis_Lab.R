# Load GEOquery library
BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("pd.hg.u133.plus.2")
BiocManager::install("hgu133plus2.db")
BiocManager::install("genefilter")
BiocManager::install("limma")
BiocManager::install("pheatmap")

library(GEOquery)

# Download pre-computed expression table
gse <- getGEO("GSE55945",GSEMatrix = T)

# Extract expression table
class(gse[[1]])
data.matrix.GEO <- exprs(gse[[1]])
boxplot(log2(data.matrix.GEO))

# Download CEL raw files
#options(timeout = max(300, getOption("timeout")))
#options(download.file.method.GEOquery = "wget")
file_paths <- getGEOSuppFiles("GSE55945",baseDir = "/home/aromanel/")

# Set working directory
setwd("/home/aromanel/GSE55945/")

# Extract files
system(paste("tar -xvf ",rownames(file_paths)[1],sep=""))

# Library for microarray data managment
library(oligo)
library("pd.hg.u133.plus.2") # Platform Design Info for Affymetrix HG-U133_Plus_2

# import CEL files containing raw probe-level data into an R AffyBatch object
celpath <- "/home/aromanel/GSE55945/"
list <- list.files(celpath,full.names=TRUE,pattern = "gz")
data <- read.celfiles(list)

list <- list[-grep("GSM1348937_110807_HGU133_PLUS_2.0_MS_36A6.CEL.gz",list)]
data <- read.celfiles(list)

list <- list[-grep("GSM1348948_011508_HGU133_PLUS_2.0_MS_36D2.CEL.gz",list)]
data <- read.celfiles(list)

# Normalize data using RMA method
data.rma <- rma(data)

# Extract annotations (if any)
pData(data.rma)
fData(data.rma)

# See microarray image
image(data[,1],main=pData(data)$index[1])

# Signal intensity profile check
hist(data[,1],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=pData(data)$index[1])
hist(data[,2],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=pData(data)$index[2])
hist(data[,3],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=pData(data)$index[3])

# Boxplot before normalization
boxplot(data,which="pm",names=pData(data)$index) 
# Boxplot after normalization
boxplot(data.rma,names=pData(data)$index)

# PCA plot
data.matrix <- exprs(data.rma)
color <- c(rep("red",12),rep("blue",7))
data.PC <- prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1:2],xlim=c(-200,200),ylim=c(-200,200),col=color,pch=19)
names <- sapply(rownames(pData(data)),function(x) strsplit(x,"\\_")[[1]][1])
text(data.PC$x[,1],data.PC$x[,2]+10,names,cex=0.7)

# Summarization by best gene probe
library(hgu133plus2.db) # Affymetrix HG-U133_Plus_2 Array annotation data
library(genefilter)

iqrs <- apply(exprs(data.rma), 1, IQR)
# select probes with largest variability
prbs <- findLargest(featureNames(data.rma), testStat = iqrs, data = "hgu133plus2.db")
dim(data.rma)
data.rma <- data.rma[prbs, ]
dim(data.rma)
ann <- AnnotationDbi::select(hgu133plus2.db, keys = featureNames(data.rma), columns = c("ENTREZID", "SYMBOL", "GENENAME"))
rownames(ann) <- ann[,1]
fData(data.rma) <- ann
data.matrix <- exprs(data.rma)

# PCA plot
color <- c(rep("red",12),rep("blue",7))
data.PC <- prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1:2],xlim=c(-200,200),ylim=c(-200,200),col=color,pch=19)

# Differential gene expression analysis (basics analysis)
genes <- fData(data.rma)$SYMBOL
annt <- c(rep("Tumor",12),rep("Normal",7))
pvals <- sapply(genes,function(x) 
   wilcox.test(data.matrix[which(genes==x),which(annt=="Tumor")],
               data.matrix[which(genes==x),which(annt=="Normal")])$p.value)
pvals.adj <- p.adjust(pvals,method = "fdr")
diff.genes <- genes[which(pvals.adj<0.05)]

# Differential gene expression analysis (standard analysis)
library(limma)
design <- model.matrix(~0+annt)
colnames(design) <- c("Normal","Tumor")
contr <- makeContrasts(Tumor-Normal,levels=design)
fit <- lmFit(data.matrix,design)
fit.contr <- eBayes(contrasts.fit(fit,contr))
etab <- topTable(fit.contr, 
                 number = 1000, 
                 adjust = "BH", 
                 p.value = 0.05, 
                 coef = "Tumor - Normal",genelist=genes)
head(etab)
diff.genes.2 <- etab$ID

data.rma.DEGs <- data.rma[which(fData(data.rma)$SYMBOL%in%diff.genes.2),]
data.matrix.DEGs <- exprs(data.rma.DEGs)

length(diff.genes)
length(diff.genes.2)
length(intersect(diff.genes,diff.genes.2))

# filter also lfc
etab <- topTable(fit.contr, number = 1000, adjust = "BH", p.value = 0.05, lfc=1, coef = "Tumor - Normal",genelist=genes)
head(etab)

# Alternative way (when only two groups)
design <- model.matrix(~annt)
colnames(design) <- c("Control","TumorVSControl")
fit <- lmFit(data.matrix,design)
fit.contr <- eBayes(fit)
etab <- topTable(fit.contr, number = 100, adjust = "BH", p.value = 0.05, lfc = 1, coef = "TumorVSControl",genelist=genes)
head(etab)

# Heatmap

dists <- as.matrix(dist(t(data.matrix.DEGs), method = "euclidean"))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))

annotation_for_heatmap <- data.frame(type = annt)
row.names(annotation_for_heatmap) <- row.names(pData(data.rma.DEGs))
ann_colors <- list(type = c(Tumor = "chartreuse4", Normal = "burlywood3"))

library(pheatmap)
pheatmap(dists, width = 90, col = (hmcol),fontsize_row = 4,annotation_legend = TRUE,
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = FALSE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small", "large")),
         main = "Prostate cancer dataset")

cols <- c(rep("chartreuse4",12),rep("burlywood3",7)) 
pal <- c("blue","white","red") 
pal <- colorRampPalette(pal)(100)
heatmap(data.matrix.DEGs,ColSideColors = cols,cexCol = 0.5,margins = c(12,4),col=pal,cexRow = 0.2)


## TASK
## Repeat all analysis with the dataset GSE3744
