#----------------------------------------------------------------#
#------------- Load required packages and functions -------------#
#----------------------------------------------------------------#

#---------------- PACKAGES ---------------#
library(Matrix)
library(Seurat)
library(sctransform)
library(tidyverse)
library(ggplot2)
library(BiocGenerics)
library(Matrix.utils)
library(reshape2)
library(cowplot)
library(ggpubr)
library(ggrastr)
library(RColorBrewer)
library(rlist)
library(pheatmap)
library(plyr)
library(dplyr)
library(openxlsx)
library(viridis)  
library(magrittr)
library(matrixStats)
library(data.table)
library(grid)
library(gtable)
library(gridExtra)
library(devtools)
library(ggvenn)
library(enrichR)
library(rcartocolor)

#--------------- FUNCTIONS ---------------#

#Smooth Norm heatmaps
flexible_normalization <-function(data_in, by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  #### if by column
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  return(output)
}

#GG color Hue
gg_color_hue <-function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Reverse %in%
'%ni%' <- Negate('%in%')


#----------------------------------------------------------------#
#------------- Quality Control + HTO demultiplexing -------------#
#----------------------------------------------------------------#

## Load RNA+HTO Lib1 raw count matrix
data <- Read10X("~/RNA_HTO_Lib1_processed/",strip.suffix = TRUE)
colnames(data[["Gene Expression"]]) <- paste0("Lib1.",colnames(data[["Gene Expression"]]))
colnames(data[["Antibody Capture"]])<- paste0("Lib1.",colnames(data[["Antibody Capture"]]))
hashing <- CreateSeuratObject(data$`Gene Expression`)
hashing[["HTO"]] <- CreateAssayObject(data$`Antibody Capture`)

hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="Seurat")

mito.genes <- grep(pattern = "^MT-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts


## Load RNA+HTO Lib2 raw count matrix
data <- Read10X("~/RNA_HTO_Lib2_processed/",strip.suffix = TRUE)
colnames(data[["Gene Expression"]]) <- paste0("Lib2.",colnames(data[["Gene Expression"]]))
colnames(data[["Antibody Capture"]])<- paste0("Lib2.",colnames(data[["Antibody Capture"]]))
hashing_2 <- CreateSeuratObject(data$`Gene Expression`)
hashing_2[["HTO"]] <- CreateAssayObject(data$`Antibody Capture`)

hashing_2 <- NormalizeData(object = hashing_2, assay = "HTO", normalization.method = "CLR")
hashing_2 <- HTODemux(hashing_2)
hashing_2 <- SetIdent(hashing_2, value="Seurat")

mito.genes <- grep(pattern = "^MT-", x = rownames(hashing_2@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(hashing_2@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing_2@assays$RNA@data == 0)/nrow(hashing_2@assays$RNA)
percent.mito <- Matrix::colSums(hashing_2@assays$RNA[mito.genes, ])/Matrix::colSums(hashing_2@assays$RNA)
percent.ribo <- Matrix::colSums(hashing_2@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing_2@assays$RNA)
hashing_2[['percent.mito']] <- percent.mito
hashing_2[['percent.ribo']] <- percent.ribo
hashing_2[['dropouts']] <- dropouts

HTO.order <- c("Hypoxie-3H_Hypoxie-6H","Hypoxie-3H_Normoxie","Hypoxie-6H_Normoxie","Hypoxie-24h_Hypoxie-3H","Hypoxie-24h_Hypoxie-6H","Hypoxie-24h_Normoxie","Normoxie","Hypoxie-3H","Hypoxie-6H","Hypoxie-24h","Negative")
HTO.colors <- c("Hypoxie-3H_Hypoxie-6H"="#F0027F","Hypoxie-3H_Normoxie"="#FDC086","Hypoxie-6H_Normoxie"="#FFFF99","Hypoxie-24h_Hypoxie-3H"="#7FC97F","Hypoxie-24h_Hypoxie-6H"="#386CB0" ,"Hypoxie-24h_Normoxie"="#BEAED4","Normoxie"="#3a8bc9","Hypoxie-3H"="#f2c0b8","Hypoxie-6H"="#de5b5b","Hypoxie-24h"="#7b3735","Negative"="#c4c3c2")
hashing$HTO_classification <- factor(hashing$HTO_classification,levels = HTO.order)
hashing_2$HTO_classification <- factor(hashing_2$HTO_classification,levels = HTO.order)
VlnPlot(hashing,group.by = "HTO_classification" ,features = c( "nCount_RNA"),y.max = 60000,pt.size = 0,cols = HTO.colors)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3))
VlnPlot(hashing_2,group.by = "HTO_classification" ,features = c( "nCount_RNA"),y.max = 60000,pt.size = 0,cols = HTO.colors)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3))

hashing$library <- "1"
hashing_2$library <- "2"
# write.table(gsub(".*\\.","",colnames(hashing)),"barcodes_lib1.tsv",row.names = F,col.names = F,quote = F)
# write.table(gsub(".*\\.","",colnames(hashing_2)),"barcodes_lib2.tsv",row.names = F,col.names = F,quote = F)


## Merging both library
all(rownames(hashing)==rownames(hashing_2))
hashing <- merge(hashing,hashing_2)

#HTO's repartition per library
metadata <- hashing@meta.data 
metadata$HTO_classification <- factor(metadata$HTO_classification, levels = HTO.order)
metadata = as.data.frame.matrix(table(metadata$library, metadata$HTO_classification))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$library = rownames(metadata)
metadata = gather(metadata, HTO_classification, percentage, "Hypoxie-3H_Hypoxie-6H":"Negative")
metadata$HTO_classification <- factor(metadata$HTO_classification, levels =HTO.order)

ggplot(metadata, aes(x = library, y = percentage, fill = HTO_classification)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = HTO.colors) +
  ggtitle("Myeloid cells relative proportions")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 0),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size =14),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



##Processing gRNA librairies with CITE-seq Count
# CITE-seq-Count -R1 ~/FASTQ/Crop_1_S3_L001_R1_001.fastq.gz -R2 ~/FASTQ/Crop_1_S3_L001_R2_001.fastq.gz -t ~/gRNA_sequences_CROPseq2.csv -cbf 1 -cbl 16 -umif 17 -umil 28 -trim 23 -cells 10000 -wl "/home/truchi/CRISPRi_2/barcodes_lib1.tsv" -o ~/gRNA_counts/Lib1_L001/
# CITE-seq-Count -R1 ~/FASTQ/Crop_1_S3_L002_R1_001.fastq.gz -R2 ~/FASTQ/Crop_1_S3_L002_R2_001.fastq.gz -t ~/gRNA_sequences_CROPseq2.csv -cbf 1 -cbl 16 -umif 17 -umil 28 -trim 23 -cells 10000 -wl "/home/truchi/CRISPRi_2/barcodes_lib1.tsv" -o ~/gRNA_counts/Lib1_L002/
# CITE-seq-Count -R1 ~/FASTQ/Crop_2_S4_L001_R1_001.fastq.gz -R2 ~/FASTQ/Crop_2_S4_L001_R2_001.fastq.gz -t ~/gRNA_sequences_CROPseq2.csv -cbf 1 -cbl 16 -umif 17 -umil 28 -trim 23 -cells 10000 -wl "/home/truchi/CRISPRi_2/barcodes_lib2.tsv" -o ~/gRNA_counts/Lib2_L001/
# CITE-seq-Count -R1 ~/FASTQ/Crop_2_S4_L002_R1_001.fastq.gz -R2 ~/FASTQ/Crop_2_S4_L002_R2_001.fastq.gz -t ~/gRNA_sequences_CROPseq2.csv -cbf 1 -cbl 16 -umif 17 -umil 28 -trim 23 -cells 10000 -wl "/home/truchi/CRISPRi_2/barcodes_lib2.tsv" -o ~/gRNA_counts/Lib2_L002/


### Load gRNA raw count matrix (4 parts)
gRNA.11 <- Read10X("~/gRNA_Lib1_L001_processed/", gene.column=1)
colnames(gRNA.11) <- paste0("Lib1.",colnames(gRNA.11))
rownames(gRNA.11) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", rownames(gRNA.11))
lib1.1 <- data.frame(lib1.1=Matrix::rowSums(gRNA.11)) %>% mutate(gRNA=rownames(gRNA.11))
data.frame(Sum=Matrix::rowSums(gRNA.11))["unmapped",]/colSums(data.frame(Sum=Matrix::rowSums(gRNA.11)))*100
#29.5% of the UMIs are unmapped

gRNA.12 <- Read10X("~/gRNA_Lib1_L002_processed/", gene.column=1)
colnames(gRNA.12) <- paste0("Lib1.",colnames(gRNA.12))
rownames(gRNA.12) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", rownames(gRNA.12))
lib1.2 <- data.frame(lib1.2=Matrix::rowSums(gRNA.12)) %>% mutate(gRNA=rownames(gRNA.12))
data.frame(Sum=Matrix::rowSums(gRNA.12))["unmapped",]/colSums(data.frame(Sum=Matrix::rowSums(gRNA.12)))*100
#28.8% of the UMIs are unmapped

gRNA.21 <- Read10X("~/gRNA_Lib2_L001_processed", gene.column=1)
colnames(gRNA.21) <- paste0("Lib2.",colnames(gRNA.21))
rownames(gRNA.21) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", rownames(gRNA.21))
lib2.1 <- data.frame(lib2.1=Matrix::rowSums(gRNA.21)) %>% mutate(gRNA=rownames(gRNA.21))
data.frame(Sum=Matrix::rowSums(gRNA.21))["unmapped",]/colSums(data.frame(Sum=Matrix::rowSums(gRNA.21)))*100
#27.7% of the UMIs are unmapped

gRNA.22 <- Read10X("~/gRNA_Lib2_L002_processed/", gene.column=1)
colnames(gRNA.22) <- paste0("Lib2.",colnames(gRNA.22))
rownames(gRNA.22) <- gsub("^([^-]*-[^-]*)-.*$", "\\1", rownames(gRNA.22))
lib2.2 <- data.frame(lib2.2=Matrix::rowSums(gRNA.22)) %>% mutate(gRNA=rownames(gRNA.22))
data.frame(Sum=Matrix::rowSums(gRNA.22))["unmapped",]/colSums(data.frame(Sum=Matrix::rowSums(gRNA.22)))*100
#27.2% of the UMIs are unmapped

summary <- Reduce(merge,list(lib1.1,lib1.2,lib2.1,lib2.2))
summary <-data.frame(row.names = summary[,1],round(sweep(summary[,-1],2,colSums(summary[,-1]),`/`)*100,digits=1))


all(rownames(gRNA.11)==rownames(gRNA.12) & rownames(gRNA.11)==rownames(gRNA.22) & rownames(gRNA.21)==rownames(gRNA.22))
gRNA <- cbind(gRNA.11,gRNA.12,gRNA.21,gRNA.22)

#check if all barcodes from gRNAs are in mRNAs
length(which(colnames(gRNA)%in%colnames(hashing)))/length(colnames(gRNA))*100

#remove unmapped counts and match colnames of gRNAs and mRNAs
gRNA <- gRNA[sort(rownames(gRNA)[-which(rownames(gRNA)%in%"unmapped")]),match(colnames(hashing),colnames(gRNA))]
dim(gRNA)
gRNA[,1:5]
hashing[["gRNA"]] <- CreateAssayObject(gRNA)

#Find cells with 0 count for all guides and remove them
cells_to_remove <- data.frame(cells=colnames(gRNA),Sum=Matrix::colSums(gRNA)) %>% filter(Sum == 0) %>% pull(cells) 
length(cells_to_remove)/dim(gRNA)[2]*100
hashing$gRNA_counts <- ">0"
Idents(hashing) <- "gRNA_counts"
hashing <- SetIdent(object = hashing, cells = cells_to_remove, value = '0')
VlnPlot(hashing,features = "nCount_RNA",pt.size = 0,y.max = 80000)
hashing <- subset(hashing,cells = cells_to_remove,invert=TRUE)
dim(hashing)

hashing <- NormalizeData(object = hashing, assay = "gRNA", normalization.method = "CLR",margin = 2)

test <- t(data.frame(GetAssayData(hashing,assay = "gRNA")))
VlnPlot(hashing,assay = "gRNA" , features = rownames(hashing@assays[["gRNA"]][1:7]),group.by = "orig.ident", ncol=7, cols = "lightsteelblue3",pt.size = 0.01,same.y.lims=TRUE)
VlnPlot(hashing,assay = "gRNA" , features = rownames(hashing@assays[["gRNA"]][8:14]), ncol=5,group.by = "orig.ident", cols = "lightsteelblue3",pt.size = 0.01,same.y.lims=TRUE)


#Assign a gRNA to each cell (+ doublets and negative)
gRNA.colors <- c("Doublet"="#470606", "HIF1A-sg1"="#b53333", "HIF1A-sg2" ="#e36262","HIF2-sg5"="#eb8831","LINC00152-sg3"="#b88cb5","LUCAT1-sg3" ="#85132f","LUCAT1-sg5"="#bd4462",  
                 "MALAT1-sg1"="#1dbf9f","NEAT1-sg2"="#c9bf24","NEAT1-sg6"="#dbd470","Neg-sg1"="#c3ccf7","Neg-sg2"="#c3def7","SNHG12-sg1" ="#42a823","SNHG12-sg3"="#69bd4f","SNHG21-sg5"="#2459c9","Negative"="#d6d2d2")
gRNA.order <- c("Doublet", "HIF1A-sg1", "HIF1A-sg2" ,"HIF2-sg5","LINC00152-sg3","LUCAT1-sg3" ,"LUCAT1-sg5",  
                "MALAT1-sg1","NEAT1-sg2","NEAT1-sg6","Neg-sg1","Neg-sg2","SNHG12-sg1","SNHG12-sg3","SNHG21-sg5","Negative")

?MULTIseqDemux()
hashing <- MULTIseqDemux(hashing,assay = "gRNA",autoThresh = TRUE)
hashing$gRNA_classification <- factor(Idents(hashing),levels = gRNA.order)
Idents(hashing) <- "gRNA_classification"
hashing$HTO_classification <- factor(hashing$HTO_classification,levels = c(HTO.order))
a <- as.data.frame.matrix(table(Idents(hashing),hashing$HTO_classification))
a <- round(sweep(a,2,colSums(a),`/`)*100,digits = 2)
a <- a[match(gRNA.order,rownames(a)),match(HTO.order,colnames(a))]
a[a>10] <- 10
pheatmap(a,cluster_rows = F,cluster_cols = F,angle_col = 90,fontsize = 8)

#Check repartition for both librairies
Idents(hashing) <- "library"
sub1 <- subset(hashing,idents="1")
Idents(sub1) <- "gRNA_classification"
a <- as.data.frame.matrix(table(Idents(sub1),sub1$HTO_classification))
a <- round(sweep(a,2,colSums(a),`/`)*100,digits = 2)
a <- a[match(gRNA.order,rownames(a)),match(HTO.order,colnames(a))]
a[a>10] <- 10
pheatmap(a,cluster_rows = F,cluster_cols = F,angle_col = 90,fontsize = 8)
sub2 <- subset(hashing,idents="2")
Idents(sub2) <- "gRNA_classification"
a <- as.data.frame.matrix(table(Idents(sub2),sub2$HTO_classification))
a <- round(sweep(a,2,colSums(a),`/`)*100,digits = 2)
a <- a[match(gRNA.order,rownames(a)),match(HTO.order,colnames(a))]
a[a>10] <- 10
pheatmap(a,cluster_rows = F,cluster_cols = F,angle_col = 90,fontsize = 8)


#Select HTO's Negative cells for further assignation method
Idents(hashing) <- "HTO_classification.global"
Negative <- subset(hashing,idents="Negative")
# saveRDS(Negative,"CRISPRi2_HTOsNeg.rds")

#Keep HTO's Singlet
hashing <- subset(hashing,idents="Singlet")
dim("Singlet")
#HTO's repartition per library
metadata <- hashing@meta.data 
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels = gRNA.order)
metadata = as.data.frame.matrix(table(metadata$library, metadata$gRNA_classification))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$library = rownames(metadata)
metadata = gather(metadata, gRNA_classification, percentage, "Doublet":"Negative")
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels =gRNA.order)

ggplot(metadata, aes(x = library, y = percentage, fill = gRNA_classification)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = gRNA.colors) +
  ggtitle("gRNA's assignation")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 0),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size =14),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

HTO.order <- c("Normoxie","Hypoxie-3H","Hypoxie-6H","Hypoxie-24h")
metadata <- hashing@meta.data 
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels = gRNA.order)
metadata$HTO_classification <- factor(metadata$HTO_classification,levels = HTO.order)
metadata = as.data.frame.matrix(table(metadata$HTO_classification, metadata$gRNA_classification))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$HTO_classification = rownames(metadata)
metadata = gather(metadata, gRNA_classification, percentage, "Doublet":"Negative")
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels = gRNA.order)
metadata$HTO_classification <- factor(metadata$HTO_classification,levels = HTO.order)
hashing$HTO_classification <- factor(hashing$HTO_classification,levels = HTO.order)

ggplot(metadata, aes(x = HTO_classification, y = percentage, fill = gRNA_classification)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = gRNA.colors) +
  ggtitle("gRNA's assignation")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size =14),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
a <- as.data.frame.matrix(table(hashing$gRNA_classification,hashing$HTO_classification))

Idents(hashing) <- "gRNA_classification"
VlnPlot(hashing,assay = "gRNA" , features = "SNHG21-sg5",pt.size = 0,ncol = 1,cols = gRNA.colors,y.max = 5)+FontSize(x.text = 8,y.text = 8,main = 10)+NoLegend()

cell.order <- hashing@meta.data %>% mutate(barcode=colnames(hashing)) %>%  arrange(gRNA_classification,HTO_classification) 
annotation_col <- data.frame(row.names = cell.order$barcode ,condition=cell.order$HTO_classification,gRNA_classification = cell.order$gRNA_classification)

#This is the heatmap of the [Figure_2B]
pdf("HTOs_heatmap.pdf",width = 10,height = 4)
pheatmap::pheatmap(data.frame(GetAssayData(hashing,assay = "gRNA"))[,match(cell.order %>% pull(barcode),colnames(GetAssayData(hashing,assay = "gRNA")))], 
                   cluster_rows = F,legend = TRUE, cluster_cols = F,annotation_col = annotation_col,show_colnames=F,annotation_legend = 1,fontsize_row=8,annotation_colors = list(gRNA_classification=gRNA.colors,condition=HTO.colors),
                   color = colorRampPalette(inferno(5))(100))
dev.off()

View(as.data.frame.matrix(table(hashing$gRNA_classification, hashing$HTO_classification)))

DefaultAssay(hashing) <- "RNA"
VlnPlot(hashing,  features = c("nCount_RNA"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,same.y.lims=F,ncol = 1)+NoLegend() +labs(y = "nCount_RNA")+FontSize(x.text = 8,y.text = 8,main = 0)
VlnPlot(hashing,  features = c("nFeature_RNA"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,same.y.lims=F,ncol = 1)+NoLegend()+labs(y = "nFeature_RNA")+FontSize(x.text = 8,y.text = 8,main = 0)

#saveRDS(hashing,"CRISPRi2_allgRNAs.rds")

#-------------------------------------------------------------------------------------------------#
#--------------------- Remove low-UMIs cells with an assigned HTO and gRNA -----------------------#
#-------------------------------------------------------------------------------------------------#

#hashing <- readRDS("CRISPRi2_allgRNAs.rds")
hashing <- hashing %>% subset(idents=c("Doublet","Negative"),invert=T)
hashing$target <- gsub("-.*","",hashing$gRNA_classification)

hashing[["RNA_raw"]] <- hashing[["RNA"]]
DefaultAssay(hashing) <- "RNA"
hashing <- NormalizeData(object = hashing,assay = "RNA") %>% ScaleData()
targets <- c("HIF1A","EPAS1","CYTOR","LUCAT1","MALAT1","NEAT1","SNHG12","SNHG21")
VlnPlot(hashing,  features = c("nCount_RNA","nFeature_RNA","percent.mito","dropouts"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0.01,same.y.lims=F,ncol = 2)+NoLegend()
VlnPlot(hashing,  features = c("SNHG21"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,same.y.lims=F,ncol = 1)+NoLegend()+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "HIF1A expression")

metadata <- hashing@meta.data %>% dplyr::select(target,gRNA_classification) %>% unique() %>% filter(target%ni%"Neg")
gRNA.order <- c("HIF1A-sg1", "HIF1A-sg2" ,"HIF2-sg5","LINC00152-sg3","LUCAT1-sg3" ,"LUCAT1-sg5",  
                "MALAT1-sg1","NEAT1-sg2","NEAT1-sg6","SNHG12-sg1","SNHG12-sg3","SNHG21-sg5")
metadata <- metadata[match(gRNA.order,metadata$gRNA_classification),]
metadata$target_gene <- c("HIF1A","HIF1A","EPAS1","CYTOR","LUCAT1","LUCAT1","MALAT1","NEAT1","NEAT1","SNHG12","SNHG12","SNHG21")

meta <- data.frame(gRNA_classification=gRNA.order)
for (cond in 1:length(HTO.order)) {
  Idents(hashing) <- "HTO_classification"
  SUB <- subset(hashing,idents=HTO.order[cond])
  for (i in 1:nrow(metadata)) {
    Idents(SUB) <- "gRNA_classification"
    sub <- subset(SUB,idents=c(as.character(metadata$gRNA_classification[i])))
    table(sub$gRNA_classification)
    expr_1 <- mean(sub@assays$RNA[metadata$target_gene[i],])
    Idents(SUB) <- "target"
    sub2 <- subset(SUB,idents="Neg")
    table(sub2$target)
    expr_2 <- mean(sub2@assays$RNA[metadata$target_gene[i],])
    metadata$pct_inhib[i] <-( expr_1-expr_2 )/expr_2*100
    metadata <- metadata %>% dplyr::select(target_gene,gRNA_classification,pct_inhib) 
  }
  meta <- merge(meta,metadata,by="gRNA_classification")
}
meta <- meta[,c(1,3,5,7,9)]
colnames(meta) <- c("gRNA_classification",HTO.order)

sub1 <- subset(hashing,idents="Normoxie")
sub2 <- subset(hashing,idents="Hypoxie-3H")
sub3 <- subset(hashing,idents="Hypoxie-6H")
sub4 <- subset(hashing,idents="Hypoxie-24h")

df <- data.frame(row.names=HTO.order,
                 nCount_RNA=c(mean(sub1$nCount_RNA),mean(sub2$nCount_RNA),mean(sub3$nCount_RNA),mean(sub4$nCount_RNA)),
                 nFeature_RNA=c(mean(sub1$nFeature_RNA),mean(sub2$nFeature_RNA),mean(sub3$nFeature_RNA),mean(sub4$nFeature_RNA)),
                 pct.mito=c(mean(sub1$percent.mito),mean(sub2$percent.mito),mean(sub3$percent.mito),mean(sub4$percent.mito)))


ggarrange(VlnPlot(sub1,  features = c("EPAS1"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,y.max = 2.2,ncol = 1)+NoLegend()+
            geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "EPAS1 expression"),
          VlnPlot(sub2,  features = c("EPAS1"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,y.max = 2.2,ncol = 1)+NoLegend()+
            geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "EPAS1 expression"),
          VlnPlot(sub3,  features = c("EPAS1"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,y.max = 2.2,ncol = 1)+NoLegend()+
            geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "EPAS1 expression"),
          VlnPlot(sub4,  features = c("EPAS1"),group.by = "gRNA_classification",cols = gRNA.colors,pt.size = 0,y.max = 2.2,ncol = 1)+NoLegend()+
            geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "EPAS1 expression"),ncol=4,nrow=1)


HTO.colors <- c("Normoxie"="#3a8bc9","Hypoxie-3H"="#f2c0b8","Hypoxie-6H"="#de5b5b","Hypoxie-24h"="#7b3735")
target.colors <- c("HIF1A" ="#e36262","HIF2"="#eb8831","LINC00152"="#b88cb5","LUCAT1" ="#85132f", "MALAT1"="#1dbf9f","NEAT1"="#c9bf24","Neg"="#c3ccf7","SNHG12" ="#42a823","SNHG21"="#2459c9")

DefaultAssay(hashing) <- "RNA_raw"
hashing <-hashing %>% SCTransform()
hashing <- RunPCA(object = hashing,assay = "SCT")
hashing <- RunUMAP(object = hashing, dims = 1:50)
hashing <- FindNeighbors(hashing, dims = 1:50, verbose = FALSE,k.param = 10)
hashing <- FindClusters(hashing, resolution = 0.3, verbose = FALSE)
DimPlot(hashing)

DefaultAssay(hashing) <- "RNA"
a <- FindAllMarkers(hashing,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
DotPlot(hashing,features = a %>% group_by(cluster) %>% top_n(10,avg_log2FC) %>% pull(gene),cols = "RdBu",dot.scale = 4) + RotatedAxis() +FontSize(x.text = 8,y.text = 8)

metadata <- hashing@meta.data 
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels = gRNA.order)
metadata = as.data.frame.matrix(table(metadata$seurat_clusters, metadata$gRNA_classification))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$seurat_clusters = rownames(metadata)
metadata = gather(metadata, gRNA_classification, percentage, "HIF1A-sg1":"SNHG21-sg5")
metadata$gRNA_classification <- factor(metadata$gRNA_classification, levels = gRNA.order)

ggplot(metadata, aes(x = seurat_clusters, y = percentage, fill = gRNA_classification)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = gRNA.colors) +
  ggtitle("gRNA's assignation")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size =14),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

metadata <- hashing@meta.data 
metadata$HTO_classification <- factor(metadata$HTO_classification, levels = HTO.order)
metadata = as.data.frame.matrix(table(metadata$seurat_clusters, metadata$HTO_classification))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$seurat_clusters = rownames(metadata)
metadata = gather(metadata, HTO_classification, percentage, "Normoxie":"Hypoxie-24h")
metadata$HTO_classification <- factor(metadata$HTO_classification, levels = HTO.order)

ggplot(metadata, aes(x = seurat_clusters, y = percentage, fill = HTO_classification)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = HTO.colors) +
  ggtitle("gRNA's assignation")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size =14),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


ggarrange(FeaturePlot(hashing,features = c("nCount_RNA"),ncol = 1 )+scale_color_viridis(),FeaturePlot(hashing,features = c("nFeature_RNA"),ncol = 1 )+scale_color_viridis(),
          FeaturePlot(hashing,features = c("percent.mito"),ncol = 1 )+scale_color_viridis(),FeaturePlot(hashing,features = c("dropouts"),ncol = 1 )+scale_color_viridis(),nrow = 2,ncol = 2)

DimPlot(hashing, group.by = "gRNA_classification", label = F, pt.size = 0.2, reduction = "umap",cols = gRNA.colors ) + ggtitle("gRNA_classification") 
DimPlot(hashing, group.by = "HTO_classification", label = F, pt.size = 0.2, reduction = "umap",cols = HTO.colors ) + ggtitle("HTO_classification") 

VlnPlot(hashing,group.by = "seurat_clusters" ,features = c( "nCount_RNA"),y.max = 60000,pt.size = 0)
VlnPlot(hashing,group.by = "seurat_clusters" ,features = c( "percent.mito"),pt.size = 0)

#Remove lowQC cluster
hashing <- subset(hashing,idents=c(5,6),invert=TRUE)
DefaultAssay(hashing) <- "RNA_raw"
hashing <-hashing %>% SCTransform()
hashing <- RunPCA(object = hashing,assay = "SCT")
hashing <- RunUMAP(object = hashing, dims = 1:50)
hashing <- FindNeighbors(hashing, dims = 1:50, verbose = FALSE,k.param = 10)
hashing <- FindClusters(hashing, resolution = 0.3, verbose = FALSE)
DimPlot(hashing,label = T)
DimPlot(hashing, group.by = "gRNA_classification", label = F, pt.size = 0.2, reduction = "umap",cols = gRNA.colors ) + ggtitle("gRNA_classification") 
DimPlot(hashing, group.by = "HTO_classification", label = F, pt.size = 0.2, reduction = "umap",cols = HTO.colors ) + ggtitle("HTO_classification") 


plot.list <- list()
for (i in 1:length(levels(hashing$gRNA_classification))) {
  gRNA.colors.bis <- setNames(rep("gainsboro",length(levels(hashing$gRNA_classification))),as.character(levels(hashing$gRNA_classification)))
  gRNA.colors.bis[i] <- "red" 
  p <- DimPlot(hashing, group.by = "gRNA_classification", label = F, pt.size = 0.2, reduction = "umap",cols =gRNA.colors.bis ) + 
    ggtitle(levels(hashing$gRNA_classification)[i])  +NoLegend()+NoAxes()
  plot.list <- list.append(plot.list,p)
}
plot_grid(plotlist = plot.list[c(1:8)],ncol = 4)
plot_grid(plotlist = plot.list[c(9:16)],ncol = 4)

Hypoxia.genes <- c("ALDOA","ANLN","BNIP3","CA9","CDKN3","CHCHD2","CTSL2","DDIT4","ENO1","GAPDH","GPI","KIF20A","KIF4A","LDHA","MCTS1","MRPL13","MRPL15","MRPS17","NDRG1","PFKP","PGK1","PSRC1","SHCBP1","SLC16A1","SLC2A1","TPI1","TUBA1C")
hashing <- AddModuleScore(hashing,features = list(Hypoxia.genes),name = "Hypoxic_score_")
FeaturePlot(hashing,features = c("Hypoxic_score_1"),ncol = 1,cols = colorRampPalette(rev(brewer.pal(n = 8, name ="Spectral")))(6))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hashing <- CellCycleScoring(hashing, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(hashing,group.by = "Phase")+ggtitle("Cell_Cycle_Score")+theme(plot.title = element_text(hjust = 0.5)) 
FeaturePlot(hashing,features = c("HIF1A"),ncol = 1,cols = colorRampPalette(rev(brewer.pal(n = 8, name ="Spectral")))(6))

DefaultAssay(hashing ) <- "RNA"
hashing[["SCT"]] <- NULL
hashing$gRNA_classification <- factor(hashing$gRNA_classification,levels= c("HIF1A-sg1", "HIF1A-sg2","HIF2-sg5","LINC00152-sg3","LUCAT1-sg3","LUCAT1-sg5","MALAT1-sg1","NEAT1-sg2","NEAT1-sg6","SNHG12-sg1","SNHG12-sg3","SNHG21-sg5","Neg-sg1","Neg-sg2"))
# saveRDS(hashing,file = "CRISPRi2_filtered_gRNAs.rds")


gRNA.colors <- c("HIF1A-sg1"="#b53333", "HIF1A-sg2" ="#e36262","HIF2-sg5"="#eb8831","LINC00152-sg3"="#b88cb5","LUCAT1-sg3" ="#85132f","LUCAT1-sg5"="#bd4462",  
                 "MALAT1-sg1"="#1dbf9f","NEAT1-sg2"="#c9bf24","NEAT1-sg6"="#dbd470","SNHG12-sg1" ="#42a823","SNHG12-sg3"="#69bd4f","SNHG21-sg5"="#2459c9","Neg-sg1"="#c3ccf7","Neg-sg2"="#c3def7")
gRNA.order <- c("HIF1A-sg1", "HIF1A-sg2" ,"HIF2-sg5","LINC00152-sg3","LUCAT1-sg3" ,"LUCAT1-sg5",  
                "MALAT1-sg1","NEAT1-sg2","NEAT1-sg6","SNHG12-sg1","SNHG12-sg3","SNHG21-sg5","Neg-sg1","Neg-sg2")
HTO.colors <- c("Normoxie"="#3a8bc9","Hypoxie-3H"="#f2c0b8","Hypoxie-6H"="#de5b5b","Hypoxie-24h"="#7b3735")
target.colors <- c("HIF1A" ="#e36262","HIF2"="#eb8831","LINC00152"="#b88cb5","LUCAT1" ="#85132f", "MALAT1"="#1dbf9f","NEAT1"="#c9bf24","Neg"="#c3ccf7","SNHG12" ="#42a823","SNHG21"="#2459c9")
targets <- c("HIF1A","EPAS1","CYTOR","LUCAT1","MALAT1","NEAT1","SNHG12","SNHG21")

# hashing <- readRDS(file = "CRISPRi2_filtered_gRNAs.rds")

ggarrange(VlnPlot(hashing,group.by = "HTO_classification" ,features = "HIF1A",pt.size = 0,cols = HTO.colors)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "EPAS1",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "CYTOR",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "LUCAT1",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "MALAT1",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "NEAT1",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "SNHG12",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),
          VlnPlot(hashing,group.by = "HTO_classification" ,features = "SNHG21",pt.size = 0,cols = HTO.colors,ncol = 1)+FontSize(x.text=8, y.text=9,x.title=0, y.title=0),ncol = 8,nrow = 1,common.legend = T)


hashing$group <- paste0(hashing$HTO_classification,"_",hashing$gRNA_classification)
Idents(hashing) <- "group"
DefaultAssay(hashing ) <- "RNA"
data <- AverageExpression(hashing,assays = "RNA",return.seurat = T)
data <- GetAssayData(data)[targets,]
annotation <- data.frame(id = colnames(data)) 
annotation <- annotation %>% mutate(condition=gsub("_.*","",id),gRNAs=gsub(".*_","",id))
annotation$condition <- factor(annotation$condition,levels = c("Normoxie","Hypoxie-3H","Hypoxie-6H","Hypoxie-24h"))
annotation$gRNAs <- factor(annotation$gRNAs,levels = gRNA.order)
annotation <- annotation %>% arrange(gRNAs ,condition) %>% column_to_rownames(.,"id")
data <- data[,match(rownames(annotation),colnames(data))]
data <- flexible_normalization(data)
data[data< -2] <- -2
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

#This is the heatmap of the [Figure_2C]
pdf("heatmap.pdf",width = 10,height = 3)
print(pheatmap(data,annotation = annotation,color = myColor,cluster_rows = F,gaps_col = c(4,8,12,16,20,24,28,32,36,40,44,48,52),cluster_cols = F,treeheight_col = 0,show_colnames=F,breaks = Breaks,annotation_legend = 0,annotation_colors = list("gRNAs"=gRNA.colors,"condition"=HTO.colors)))
dev.off()


#-------------------------------------------------------------------------------------------------------------#
#-------------------------- Prepare input data for autoencoder without gRNA counts ---------------------------#
#-------------------------------------------------------------------------------------------------------------#

hashing <- readRDS(file = "/data/truchi_data/Rdata_objects/CRISPRi2_filtered_gRNAs.rds")
hashing$barcode <- colnames(hashing)
Idents(hashing) <- "HTO_classification"
hashing_Hx3h <- subset(hashing,idents="Hypoxie-3H")
hashing_Hx6h <- subset(hashing,idents="Hypoxie-6H")
hashing_Hx24h <- subset(hashing,idents="Hypoxie-24h")
hashing_Nx <- subset(hashing,idents="Normoxie") 
gRNA.colors <- c("HIF1A-sg1"="#b53333", "HIF1A-sg2" ="#e36262","HIF2-sg5"="#eb8831","LINC00152-sg3"="#b88cb5","LUCAT1-sg3" ="#85132f","LUCAT1-sg5"="#bd4462",  
                 "MALAT1-sg1"="#1dbf9f","NEAT1-sg2"="#c9bf24","NEAT1-sg6"="#dbd470","Neg-sg1"="#c3ccf7","Neg-sg2"="#c3def7","SNHG12-sg1" ="#42a823","SNHG12-sg3"="#69bd4f","SNHG21-sg5"="#2459c9")
target.colors <- c("HIF1A" ="#e36262","HIF2"="#eb8831","LINC00152"="#b88cb5","LUCAT1" ="#85132f", "MALAT1"="#1dbf9f","NEAT1"="#c9bf24","Neg"="#c3ccf7","SNHG12" ="#42a823","SNHG21"="#2459c9")


targets <- c("HIF1A","HIF2","LINC00152","LUCAT1","MALAT1","NEAT1","SNHG12","SNHG21")
list.gRNA <- list("HIF1A"=c("HIF1A-sg1", "HIF1A-sg2"),"HIF2"=c("HIF2-sg5"),"LINC00152"=c("LINC00152-sg3"),"LUCAT1"=c("LUCAT1-sg3" ,"LUCAT1-sg5"),"MALAT1"=c("MALAT1-sg1"),"NEAT1"=c("NEAT1-sg2","NEAT1-sg6"),"SNHG12"=c("SNHG12-sg1","SNHG12-sg3"),"SNHG21"=c("SNHG21-sg5"))
list.barcode_Nx<- list.barcode_Hx3h <- list.barcode_Hx6h <- list.barcode_Hx24h <- list()

for (i in 1:length(targets)) {
  #Normoxie
  Idents(hashing_Nx) <- "target"
  DefaultAssay(hashing_Nx) <- "RNA_raw"
  sub <- subset(hashing_Nx,idents=c("Neg", targets[i]))
  data_raw <- GetAssayData(sub)
  
  cells <- sub@meta.data %>% dplyr::select(barcode,target) %>% dplyr::rename(Label=target)  %>% as.data.table()
  cells[ , Name := 1:.N , by = c("Label") ]
  cells$Label <- factor(cells$Label,levels = c("Neg", targets[i]))
  cells<- cells[order(cells$Label),]
  list.barcode_Nx <- list.append(list.barcode_Nx,cells)
  genes <- c(names(sort(rowSums(data_raw),decreasing = TRUE)[1:10000]))
  data_raw <- rbind(data_raw)
  
  dataframe <- data_raw[match(genes,rownames(data_raw)),match(cells$barcode,colnames(data_raw))] %>% as.matrix() %>% as.data.frame()
  cells$Name <- cells$barcode
  df <- cells %>% column_to_rownames("barcode") %>% dplyr::select(Name,Label) 
  df$Label <- ifelse(df$Label=="Neg",1,2)
  df <- rbind(df%>% t() %>% as.data.frame(),dataframe) 
  print(paste("Writing",targets[i],"matrix in Normoxia",sep = " "))
  write.table(df,paste0("~/datas/Nx_",targets[i],".nogRNA.csv"),sep = ";",col.names = F,quote = F)
  
  #Hypoxie 3h
  Idents(hashing_Hx3h) <- "target"
  DefaultAssay(hashing_Hx3h) <- "RNA_raw"
  sub <- subset(hashing_Hx3h,idents=c("Neg", targets[i]))
  data_raw <- GetAssayData(sub)
  
  cells <- sub@meta.data %>% dplyr::select(barcode,target) %>% dplyr::rename(Label=target)  %>% as.data.table()
  cells[ , Name := 1:.N , by = c("Label") ]
  cells$Label <- factor(cells$Label,levels = c("Neg", targets[i]))
  cells<- cells[order(cells$Label),]
  list.barcode_Hx3h <- list.append(list.barcode_Hx3h,cells)
  genes <- c(names(sort(rowSums(data_raw),decreasing = TRUE)[1:10000]))
  data_raw <- rbind(data_raw)
  
  dataframe <- data_raw[match(genes,rownames(data_raw)),match(cells$barcode,colnames(data_raw))] %>% as.matrix() %>% as.data.frame()
  cells$Name <- cells$barcode
  df <- cells %>% column_to_rownames("barcode") %>% dplyr::select(Name,Label) 
  df$Label <- ifelse(df$Label=="Neg",1,2)
  df <- rbind(df%>% t() %>% as.data.frame(),dataframe) 
  print(paste("Writing",targets[i],"matrix in Hypoxia 3h",sep = " "))
  write.table(df,paste0("~/datas/Hx3h_",targets[i],".nogRNA.csv"),sep = ";",col.names = F,quote = F)
  
  #Hypoxie 6h
  Idents(hashing_Hx6h) <- "target"
  DefaultAssay(hashing_Hx6h) <- "RNA_raw"
  sub <- subset(hashing_Hx6h,idents=c("Neg", targets[i]))
  data_raw <- GetAssayData(sub)
  
  cells <- sub@meta.data %>% dplyr::select(barcode,target) %>% dplyr::rename(Label=target)  %>% as.data.table()
  cells[ , Name := 1:.N , by = c("Label") ]
  cells$Label <- factor(cells$Label,levels = c("Neg", targets[i]))
  cells<- cells[order(cells$Label),]
  list.barcode_Hx6h <- list.append(list.barcode_Hx6h,cells)
  genes <- c(names(sort(rowSums(data_raw),decreasing = TRUE)[1:10000]))
  data_raw <- rbind(data_raw)
  
  dataframe <- data_raw[match(genes,rownames(data_raw)),match(cells$barcode,colnames(data_raw))] %>% as.matrix() %>% as.data.frame()
  cells$Name <- cells$barcode
  df <- cells %>% column_to_rownames("barcode") %>% dplyr::select(Name,Label) 
  df$Label <- ifelse(df$Label=="Neg",1,2)
  df <- rbind(df%>% t() %>% as.data.frame(),dataframe) 
  print(paste("Writing",targets[i],"matrix in Hypoxia 6h",sep = " "))
  write.table(df,paste0("~/datas/Hx6h_",targets[i],".nogRNA.csv"),sep = ";",col.names = F,quote = F)
  
  #Hypoxie 24h
  Idents(hashing_Hx24h) <- "target"
  DefaultAssay(hashing_Hx24h) <- "RNA_raw"
  sub <- subset(hashing_Hx24h,idents=c("Neg", targets[i]))
  data_raw <- GetAssayData(sub)
  
  cells <- sub@meta.data %>% dplyr::select(barcode,target) %>% dplyr::rename(Label=target)  %>% as.data.table()
  cells[ , Name := 1:.N , by = c("Label") ]
  cells$Label <- factor(cells$Label,levels = c("Neg", targets[i]))
  cells<- cells[order(cells$Label),]
  list.barcode_Hx24h <- list.append(list.barcode_Hx24h,cells)
  genes <- c(names(sort(rowSums(data_raw),decreasing = TRUE)[1:10000]))
  data_raw <- rbind(data_raw)
  
  dataframe <- data_raw[match(genes,rownames(data_raw)),match(cells$barcode,colnames(data_raw))] %>% as.matrix() %>% as.data.frame()
  cells$Name <- cells$barcode
  df <- cells %>% column_to_rownames("barcode") %>% dplyr::select(Name,Label) 
  df$Label <- ifelse(df$Label=="Neg",1,2)
  df <- rbind(df%>% t() %>% as.data.frame(),dataframe) 
  print(paste("Writing",targets[i],"matrix in Hypoxia 24h",sep = " "))
  write.table(df,paste0("~/datas/Hx24h_",targets[i],".nogRNA.csv"),sep = ";",col.names = F,quote = F)
}

# Save the list of barcode to link the SSAE outputs with the seurat object
names(list.barcode_Nx) <- paste0(targets,"_Nx")
names(list.barcode_Hx3h) <- paste0(targets,"_Hx3h")
names(list.barcode_Hx6h) <- paste0(targets,"_Hx6h")
names(list.barcode_Hx24h) <- paste0(targets,"_Hx24h")
list.barcode <- Reduce(append,list(list.barcode_Nx,list.barcode_Hx3h,list.barcode_Hx6h,list.barcode_Hx24h))
saveRDS(list.barcode,file = "~/list_barcode.rds")

VlnPlot(hashing_Hx24h, features = c("PEX1"),group.by = "target",pt.size = 0,ncol = 1)+NoLegend()+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=3)) +labs(y = "EPAS1 expression")
