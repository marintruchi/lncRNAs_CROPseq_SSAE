# Paper : Detecting subtle transcriptomic perturbations induced by lncRNAs Knock-Down in single-cell CRISPRi screening using a new sparse supervised autoencoder neural network
# Code : R script to produce the figures of the manuscript
# Author : Marin Truchi

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


#-------------------------------------------------------------------------------------------------------------#
#-------------------------------------- Read SSAE's classification outputs -----------------------------------#
#-------------------------------------------------------------------------------------------------------------#

gRNA_0.6.colors <- c("HIF1A-sg1"="#b53333", "HIF1A-sg2" ="#e36262","HIF2-sg5"="#eb8831","LINC00152-sg3"="#b88cb5","LUCAT1-sg3" ="#85132f","LUCAT1-sg5"="#bd4462",  
                     "MALAT1-sg1"="#1dbf9f","NEAT1-sg2"="#c9bf24","NEAT1-sg6"="#dbd470","Neg-sg1"="#c3ccf7","Neg-sg2"="#c3def7","SNHG12-sg1" ="#42a823","SNHG12-sg3"="#69bd4f","SNHG21-sg5"="#2459c9")
target.colors <- c("HIF1A" ="#e36262","HIF2"="#eb8831","LINC00152"="#b88cb5","LUCAT1" ="#85132f", "MALAT1"="#1dbf9f","NEAT1"="#c9bf24","Neg"="#c3ccf7","SNHG12" ="#42a823","SNHG21"="#2459c9")
targets <- c("HIF1A","HIF2","LINC00152","LUCAT1","MALAT1","NEAT1","SNHG12","SNHG21")
timepoints <- c("Nx","Hx3h","Hx6h","Hx24h")

## Load Seurat object
hashing <- readRDS(file = "~/processed_seurat_obj.rds")
hashing$barcode <- colnames(hashing)
Idents(hashing) <- "HTO_classification"

## Here, for each target and each time points, we load the classification scores obtained for each of the 3 run
## The cell is classified as Perturbed only if 
list_prtb <- list()
for (time in timepoints) {
  for ( target in targets) {
    label1 <- read.csv(paste0("~/results_stat/",time,"_",target,"_1/Labelspred_softmax.csv"),sep = ";") %>% filter(Labels==1 & Proba.class.1 > 0.5) %>% pull(Name)
    label2 <- read.csv(paste0("~/results_stat/",time,"_",target,"_2/Labelspred_softmax.csv"),sep = ";") %>% filter(Labels==1 & Proba.class.1 > 0.5) %>% pull(Name)
    label3 <- read.csv(paste0("~/results_stat/",time,"_",target,"_3/Labelspred_softmax.csv"),sep = ";") %>% filter(Labels==1 & Proba.class.1 > 0.5) %>% pull(Name)
    a <- Reduce(f=intersect,list(label1,label2,label3))
    list_prtb <- list.append(list_prtb,a)
  }
}

list_prtb <- unlist(list_prtb)
hashing$autoencoder_final_class <- ifelse(hashing$barcode %in% list_prtb,"PRTB",ifelse(hashing$target=="Neg","NT","NP"))
hashing$autoencoder_final_class <- factor(hashing$autoencoder_final_class,levels = c("PRTB","NP","NT"))
DefaultAssay(hashing) <- "RNA"
Idents(hashing) <- "target"

hashing_Hx3h <- subset(hashing,subset=HTO_classification=="Hypoxie-3H")
hashing_Hx6h <- subset(hashing,subset=HTO_classification=="Hypoxie-6H")
hashing_Hx24h <- subset(hashing,subset=HTO_classification=="Hypoxie-24h")
hashing_Nx <- subset(hashing,subset=HTO_classification=="Normoxie")

View(as.data.frame.matrix(table(hashing_Nx$target,hashing_Nx$autoencoder_final_class)))
View(as.data.frame.matrix(table(hashing_Hx3h$target,hashing_Hx3h$autoencoder_final_class)))
View(as.data.frame.matrix(table(hashing_Hx6h$target,hashing_Hx6h$autoencoder_final_class)))
View(as.data.frame.matrix(table(hashing_Hx24h$target,hashing_Hx24h$autoencoder_final_class)))


# Calculate percentage of PRTB cells for all target gene classes.
df <- prop.table(table(hashing_Nx$autoencoder_final_class, hashing_Nx$gRNA_classification),2)
print(t(table(hashing_Nx$autoencoder_final_class, hashing_Nx$gRNA_classification)))
df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "PRTB"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c( "NP", "PRTB"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-sg")[[1]][2])
df2 <- df2 %>% filter(gene%ni%"Neg")
p1 <- ggplot(df2, aes(x = guide_number, y = value*100, fill= Var1)) + geom_bar(stat= "identity") +theme_classic()+
  scale_fill_manual(values = c("#ede4d1","#3a8bc9")) + ylab("% of cells") + xlab("sgRNA")+
  theme(axis.text.x = element_text(size = 14, hjust = 1), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 12), 
        strip.text = element_text(size=10, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 11, scales = "free") +
  labs(fill = "autoencoder class") +theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12))+xlab("")+NoLegend()

df <- prop.table(table(hashing_Hx3h$autoencoder_final_class, hashing_Hx3h$gRNA_classification),2)
print(t(table(hashing_Hx3h$autoencoder_final_class, hashing_Hx3h$gRNA_classification)))
df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "PRTB"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c( "NP", "PRTB"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-sg")[[1]][2])
df2 <- df2 %>% filter(gene%ni%"Neg")
p2 <- ggplot(df2, aes(x = guide_number, y = value*100, fill= Var1)) + geom_bar(stat= "identity") +theme_classic()+
  scale_fill_manual(values = c("#ede4d1","#f2c0b8")) + ylab("% of cells") + xlab("sgRNA")+
  theme(axis.text.x = element_text(size = 14, hjust = 1), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 12), 
        strip.text = element_text(size=10, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 11, scales = "free") +
  labs(fill = "autoencoder class") +theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12))+xlab("")+NoLegend()

df <- prop.table(table(hashing_Hx6h$autoencoder_final_class, hashing_Hx6h$gRNA_classification),2)
print(t(table(hashing_Hx6h$autoencoder_final_class, hashing_Hx6h$gRNA_classification)))
df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "PRTB"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c( "NP", "PRTB"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-sg")[[1]][2])
df2 <- df2 %>% filter(gene%ni%"Neg")
p3 <- ggplot(df2, aes(x = guide_number, y = value*100, fill= Var1)) + geom_bar(stat= "identity") +theme_classic()+
  scale_fill_manual(values = c("#ede4d1","#de5b5b")) + ylab("% of cells") + xlab("sgRNA") + 
  theme(axis.text.x = element_text(size = 14, hjust = 1), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 12), 
        strip.text = element_text(size=10, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 11, scales = "free") +
  labs(fill = "autoencoder class") +theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12))+xlab("")+NoLegend()

df <- prop.table(table(hashing_Hx24h$autoencoder_final_class, hashing_Hx24h$gRNA_classification),2)
print(t(table(hashing_Hx24h$autoencoder_final_class, hashing_Hx24h$gRNA_classification)))
df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "PRTB"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c( "NP", "PRTB"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "-sg")[[1]][2])
df2 <- df2 %>% filter(gene%ni%"Neg")
p4 <- ggplot(df2, aes(x = guide_number, y = value*100, fill= Var1)) + geom_bar(stat= "identity") +theme_classic()+
  scale_fill_manual(values = c("#ede4d1","#7b3735")) + ylab("% of cells") + xlab("sgRNA") + theme(axis.text.x = element_text(size = 14, hjust = 1), 
                                                                                                  axis.text.y = element_text(size = 14), 
                                                                                                  axis.title = element_text(size = 12), 
                                                                                                  strip.text = element_text(size=10, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 11, scales = "free") +
  labs(fill = "autoencoder class") +theme(legend.title = element_text(size = 14),legend.text = element_text(size = 12))+xlab("")+NoLegend()


#This is the plot of the [Figure_4]
ggarrange(p1,p2,p3,p4,ncol = 1,nrow = 4)
ggsave("Fig_percent_PRTB.pdf", width = 25, height = 22, units = c("cm"), dpi = 200)



## Load the lists of the most discriminant features between cells targeted by a particular gRNA and control cells
## Here, for each target and each time points, we load 3 different lists from 3 independent SSAE run on the same input data (with different initialization seed for each run),
## in order to compute a mean and a standard deviation of the obtained ranks 
for (time in timepoints) {
  for ( target in targets) {
    a <- read.csv(paste0("~/results_stat/",time,"_",target,"_1/proj_l11ball_topGenes_Mean_Captum_dl_300.csv"),sep = ";") %>% mutate(rank1=row_number(-Mean))%>% dplyr::select(Features,rank1)
    b <- read.csv(paste0("~/results_stat/",time,"_",target,"_2/proj_l11ball_topGenes_Mean_Captum_dl_300.csv"),sep = ";") %>% mutate(rank2=row_number(-Mean))%>% dplyr::select(Features,rank2)
    c <- read.csv(paste0("~/results_stat/",time,"_",target,"_3/proj_l11ball_topGenes_Mean_Captum_dl_300.csv"),sep = ";") %>% mutate(rank3=row_number(-Mean))%>% dplyr::select(Features,rank3)
    final <- Reduce(inner_join,list(a,b,c)) %>% mutate(mean=rowMeans(.[,-1])) %>% arrange(mean)
    write.xlsx(final,paste0("~/results_compiled/",time,"_",target,"_compiled_ranks.xlsx"))
  }
}


## Take the top 20 best ranked genes for each target at each time point
## These are the plots of the [Figure_5A,B],[Figure_6A,B,C],[Figure_7A,B,C]
for ( target in targets) {
  a <- read.xlsx(paste0("~/results_compiled/Nx_",target,"_compiled_ranks.xlsx")) %>% slice_min(order_by = mean,n=20) %>% mutate(sd=rowSds(.[,c("rank1","rank2","rank3")] %>% as.matrix())) 
  sub1 <- subset(hashing_Nx,cells = hashing_Nx@meta.data %>% filter(target==target & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode))
  df <- FindMarkers(sub1,ident.1 = target,ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=a$Features) %>% mutate(sign=ifelse(avg_log2FC>0,"up","down")) %>% mutate(Features=rownames(.)) %>% dplyr::select(Features,sign)
  df1 <-  inner_join(a,df)%>% arrange(mean)
  df1$Features <- factor(df1$Features,levels = rev(df1$Features))
  p1 <- ggplot(df1, aes(x=Features, y=mean)) + 
    geom_bar(stat="identity", position=position_dodge(), fill="#3a8bc9", width=.8) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
    theme_minimal()+xlab("")+ylab("")+coord_flip(ylim = c(1, 100))+rotate_x_text(45)+theme(axis.text.y = element_text(hjust = 1, colour = ifelse(rev(df1$sign) == "up", "#a3070b", "#034a96")))
  p1 
  
  a <- read.xlsx(paste0("~/results_compiled/Hx3h_",target,"_compiled_ranks.xlsx")) %>% slice_min(order_by = mean,n=20) %>% mutate(sd=rowSds(.[,c("rank1","rank2","rank3")] %>% as.matrix())) %>% arrange(mean)
  sub1 <- subset(hashing_Hx3h,cells = hashing_Hx3h@meta.data %>% filter(target==target & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode))
  df <- FindMarkers(sub1,ident.1 = target,ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=a$Features) %>% mutate(sign=ifelse(avg_log2FC>0,"up","down")) %>% mutate(Features=rownames(.)) %>% dplyr::select(Features,sign)
  df2 <-  inner_join(a,df)%>% arrange(mean)
  df2$Features <- factor(df2$Features,levels = rev(df2$Features))
  p2 <- ggplot(df2, aes(x=Features, y=mean)) + 
    geom_bar(stat="identity", position=position_dodge(), fill="#f2c0b8", width=.8) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
    theme_minimal()+xlab("")+ylab("")+coord_flip(ylim = c(1, 100))+rotate_x_text(45)+theme(axis.text.y = element_text(hjust = 1, colour = ifelse(rev(df2$sign) == "up", "#a3070b", "#034a96")))
  p2 
  
  a <- read.xlsx(paste0("~/results_compiled/Hx6h_",target,"_compiled_ranks.xlsx")) %>% slice_min(order_by = mean,n=20) %>% mutate(sd=rowSds(.[,c("rank1","rank2","rank3")] %>% as.matrix())) %>% arrange(mean)
  sub1 <- subset(hashing_Hx6h,cells = hashing_Hx6h@meta.data %>% filter(target==target & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode))
  df <- FindMarkers(sub1,ident.1 = target,ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=a$Features) %>% mutate(sign=ifelse(avg_log2FC>0,"up","down")) %>% mutate(Features=rownames(.)) %>% dplyr::select(Features,sign)
  df3 <-  inner_join(a,df)%>% arrange(mean)
  df3$Features <- factor(df3$Features,levels = rev(df3$Features))
  p3 <- ggplot(df3, aes(x=Features, y=mean)) + 
    geom_bar(stat="identity", position=position_dodge(), fill="#de5b5b", width=.8) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
    theme_minimal()+xlab("")+ylab("")+coord_flip(ylim = c(1, 100))+rotate_x_text(45)+theme(axis.text.y = element_text(hjust = 1, colour = ifelse(rev(df3$sign) == "up", "#a3070b", "#034a96")))
  p3
  
  a <- read.xlsx(paste0("~/results_compiled/Hx24h_",target,"_compiled_ranks.xlsx")) %>% slice_min(order_by = mean,n=20) %>% mutate(sd=rowSds(.[,c("rank1","rank2","rank3")] %>% as.matrix())) %>% arrange(mean)
  sub1 <- subset(hashing_Hx24h,cells = hashing_Hx24h@meta.data %>% filter(target==target & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode))
  df <- FindMarkers(sub1,ident.1 = target,ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=a$Features) %>% mutate(sign=ifelse(avg_log2FC>0,"up","down")) %>% mutate(Features=rownames(.)) %>% dplyr::select(Features,sign)
  df4 <-  inner_join(a,df)%>% arrange(mean)
  df4$Features <- factor(df4$Features,levels = rev(df4$Features))
  p4 <- ggplot(df4, aes(x=Features, y=mean)) + 
    geom_bar(stat="identity", position=position_dodge(), fill="#7b3735", width=.8) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,position=position_dodge(.9))+
    theme_minimal()+xlab("")+ylab("")+coord_flip(ylim = c(1, 100))+rotate_x_text(45)+theme(axis.text.y = element_text(hjust = 1, colour = ifelse(rev(df4$sign) == "up", "#a3070b", "#034a96")))
  p4
  
  ggarrange(p1,p2,p3,p4,ncol = 4)
  ggsave(paste0("Fig_ranks_",target,".pdf"), width = 28, height = 10, units = c("cm"), dpi = 200)
  
}


## Comparison of the perturbation signatures induced by the Knock-Down of HIF1A or HIF2 at each time point
HTO.colors <- c("Normoxie"="#3a8bc9","Hypoxie-3H"="#f2c0b8","Hypoxie-6H"="#de5b5b","Hypoxie-24h"="#7b3735")

dfNx <- read.xlsx("~/results_compiled/Nx_HIF1A_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx3h <- read.xlsx("~/results_compiled/Hx3h_HIF1A_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx6h <- read.xlsx("~/results_compiled/Hx6h_HIF1A_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx24h <- read.xlsx("~/results_compiled/Hx24h_HIF1A_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
HIF1A <- c(dfNx,dfHx3h,dfHx6h,dfHx24h) %>% unique()


dfNx <- read.xlsx("~/results_compiled/Nx_HIF2_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx3h <- read.xlsx("~/results_compiled/Hx3h_HIF2_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx6h <- read.xlsx("~/results_compiled/Hx6h_HIF2_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
dfHx24h <- read.xlsx("~/results_compiled/Hx24h_HIF2_compiled_ranks.xlsx") %>% slice_min(order_by = mean,n=20) %>% pull(Features)
HIF2 <- c(dfNx,dfHx3h,dfHx6h,dfHx24h) %>% unique()

targets <- c("HIF1A","HIF2")
dea <- list()
for (i in 1:length(targets)) {
  #select only perturbed cells
  cells <- hashing_Nx@meta.data %>% filter(target==targets[i] & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode)
  sub1 <- subset(hashing_Nx,cells = cells)
  table(sub1$target)
  #perform diff. expr. analysis between PRTB and Control cells
  print(paste("Performing DEA for",targets[i],"in Normoxia",sep = " "))
  DefaultAssay(sub1) <- "RNA"
  Idents(sub1) <- "target"
  df <- FindMarkers(sub1,ident.1 = targets[i],ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=unique(c(HIF1A,HIF2)))
  df1 <- data.frame(gene=rownames(df),Nx=df$avg_log2FC) %>% dplyr::select(gene,Nx)
  
  cells <- hashing_Hx3h@meta.data %>% filter(target==targets[i] & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode)
  sub1 <- subset(hashing_Hx3h,cells = cells)
  table(sub1$target)
  #perform diff. expr. analysis between PRTB and Control cells
  print(paste("Performing DEA for",targets[i],"in Hx3h",sep = " "))
  DefaultAssay(sub1) <- "RNA"
  Idents(sub1) <- "target"
  df <- FindMarkers(sub1,ident.1 = targets[i],ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=unique(c(HIF1A,HIF2)))
  df2 <- data.frame(gene=rownames(df),Hx3h=df$avg_log2FC) %>% dplyr::select(gene,Hx3h)
  
  cells <- hashing_Hx6h@meta.data %>% filter(target==targets[i] & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode)
  sub1 <- subset(hashing_Hx6h,cells = cells)
  table(sub1$target)
  #perform diff. expr. analysis between PRTB and Control cells
  print(paste("Performing DEA for",targets[i],"in Hx6h",sep = " "))
  DefaultAssay(sub1) <- "RNA"
  Idents(sub1) <- "target"
  df <- FindMarkers(sub1,ident.1 = targets[i],ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=unique(c(HIF1A,HIF2)))
  df3 <- data.frame(gene=rownames(df),Hx6h=df$avg_log2FC)  %>% dplyr::select(gene,Hx6h)
  
  cells <- hashing_Hx24h@meta.data %>% filter(target==targets[i] & autoencoder_final_class =="PRTB" | target=="Neg") %>% pull(barcode)
  sub1 <- subset(hashing_Hx24h,cells = cells)
  table(sub1$target)
  #perform diff. expr. analysis between PRTB and Control cells
  print(paste("Performing DEA for",targets[i],"in Hx24h",sep = " "))
  DefaultAssay(sub1) <- "RNA"
  Idents(sub1) <- "target"
  df <- FindMarkers(sub1,ident.1 = targets[i],ident.2 = "Neg",logfc.threshold = 0,min.pct = 0,features=unique(c(HIF1A,HIF2)))
  df4 <- data.frame(gene=rownames(df),Hx24h=df$avg_log2FC)  %>% dplyr::select(gene,Hx24h)
  df <- Reduce(inner_join,list(df1,df2,df3,df4))%>% column_to_rownames(.,"gene")
  colnames(df) <- paste0(targets[i],"_",colnames(df))
  dea <- list.append(dea,df)
}

def <- inner_join(dea[[1]] %>% mutate(gene=rownames(.)),dea[[2]] %>% mutate(gene=rownames(.)))%>% column_to_rownames(.,"gene") %>% mutate(maximum=rowMaxs(as.matrix(abs(.)))) %>% slice_max(order_by = maximum,n=40) 
def <- def[,-9]

annotation <- data.frame(id = colnames(def)) 
annotation <- annotation %>% mutate(gRNAs=gsub("_.*","",id),condition=gsub(".*_","",id))
annotation$condition <- factor(annotation$condition,levels = c("Nx","Hx3h","Hx6h","Hx24h"))
annotation$gRNAs <- factor(annotation$gRNAs,levels = c("HIF1A","HIF2"))
annotation <- annotation %>% arrange(gRNAs ,condition) %>% column_to_rownames(.,"id")
def[def< -1] <- -1
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(def), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(def)/paletteLength, max(def), length.out=floor(paletteLength/2)))

## These is the heatmap  of the [Figure_5C]
pdf("HIF1_HIF2.pdf",width = 10,height = 3)
pheatmap(t(def),breaks = Breaks,color = myColor,angle_col =45)
dev.off()


## Comparison of the 10 first selected features between SSAE (with or without cell selection),SAE and Random Forests , for HIF2 and LUCAT1 targeted cells in hypoxia 24h


library(UpSetR)

listHIF2 <- list(SSAE_cs = c("TMEM141","IGFBP3","BNIP3","FAM162A","EPAS1","PGK1","FXYD2","ATP1B1","TESC","SLC16A3"),
                 SSAE =c("TMEM141","PGK1","IGFBP3","BNIP3","EPAS1","FAM162A","FXYD2","GPI","SLC16A3","ANXA4"),
                 SAE =c("TMEM141","IGFBP3","BNIP3","FAM162A","PGK1","EPAS1","ATP1B1","LOXL2","ALDH3A1"),
                 RF = c("IGFBP3","FTL","PGK1","TMEM141","EPAS1","BNIP3","FAM162A","EPAS1","FXYD2","ATP1B1","LOXL2","ALDH3A1","ENO1"))

## These is the Venn Diagramm  of the [Figure_8A]
ggvenn(listHIF2,show_percentage = F,
       fill_color = rev(c('#39b54a','#a4d2db','#7586d9','#815aa3')),stroke_size = 0.3, set_name_size = 3,text_size =3)
ggsave("VennDiagrammHIF2.pdf", width = 6, height = 6, dpi = 800)


listLUCAT1 <- list(SSAE_cs = c("LUCAT1","ATP6AP1","KDM5C","TMEM175","PEX1","NIT1","PHF20","ZNF181","CCDC142","AGO2"),
                   SSAE =c("LUCAT1","ATP6AP1","CCDC142","NIT1","PEX1","KDM5C","PHF20","GFRA1","CMBL","ELF1"),
                   SAE =c("LUCAT1","ATP6AP1","APP","ARF6","PEX1","KLH5","PAIP2","SF3B1","KLF13","SCCPDH"),
                   RF = c("LUCAT1","BTF","G3BP1","NKAPD1","PDCL","HSF1","MFSD1","BCLAF1","UBE2R2","RPL14"))

## These is the Venn Diagramm  of the [Figure_8B]
ggvenn(listLUCAT1,show_percentage =F,
       fill_color = rev(c('#39b54a','#a4d2db','#7586d9','#815aa3')),stroke_size = 0.3, set_name_size = 3,text_size =3)
ggsave("VennDiagrammLUCAT1.pdf", width = 6, height = 6, dpi = 800)

