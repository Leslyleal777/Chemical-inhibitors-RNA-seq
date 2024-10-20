###################################### blood asxl1 bo anaylsis info #######################################################

#21 oct 21
#angela wei
#RNAseq data set: Issy RNAseq BO and control blood
#this is Deseq2,v1.32.0 on R/4.1.1
#run sessionInfo() in console to find this information
#this script will test the effects of ASXL1 BO vs controls in blood adjusting for gender

###################################### preparing data for DESeq2 #######################################################
setwd("C:/Users/angelawei/Box/_Arboleda Lab Stuff/NGS_Sequencing/RNASeq_Results/issy_ASXL1_blood_fibro/angela_blood_DEG/allSamples_adjGender")
countdata<-read.table("Issy_ASXL1_blood_featureCounts_GeneTable_final.txt",header=TRUE,row.names = 1)
library(readxl)
coldata<-data.frame(read_excel("coldata_gender.xlsx"))
rownames(coldata)<-coldata$sample
coldata$sample<-NULL
coldata$condition[coldata$condition=="ASXL1"]<-"disease"
coldata$condition[coldata$condition=="Control"]<-"control"
countdata<-as.matrix(countdata)

###################################### running DESeq2 & preliminary results #######################################################
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = countdata,colData = coldata,design = ~sex+condition)
dds<-DESeq(dds)
plotDispEsts(dds)

#make pca
rld<-rlog(dds,blind = TRUE)
#making fancy pca with labels 
#this is to plot the points on the PCA from DESeq2 with sample names on them
pcaData<-plotPCA(rld,intgroup=c("condition"),returnData=TRUE)
percentVar<-round(100*attr(pcaData,"percentVar"))
pcaData<-merge(x=pcaData,y=coldata,by=0)
library(ggplot2)
library(ggrepel)
ggplot(pcaData, aes(PC1, PC2, color=condition.x, shape=sex,label=name)) +
  geom_point(size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text_repel(aes(label = name),color="black")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#diseazse group slightly separates out on PC2
#pc1 seems to separate by gender; not adjusting for gender in this analysis
#issy already checked 1723 for mutation also is older than other

#get more than PC1 and PC2
#DESeq2 defaults to using the top 500 most variable gene 
morePCs<-function (object, intgroup = "condition", ntop = 500, returnData = FALSE,firstPC,secondPC) 
{
  #object is the log2 transformed version of the gene count table 
  rv <- rowVars(assay(object))
  #pick out the ntop (default 500) genes with the most variable expression
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  #subset the log2 transformed values to the ntop most variable gene expression and compute PCs
  pca <- prcomp(t(assay(object)[select, ]))
  #this variable is the percent explained by the chosen PCs
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, firstPC], PC2 = pca$x[, secondPC], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[firstPC:secondPC]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", 
                              color = "group")) + geom_point(size = 3) + xlab(paste0("PC",firstPC,": ", 
                                                                                     round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC",secondPC,": ", 
                                                                                                                                              round(percentVar[2] * 100), "% variance")) + coord_fixed()
}
# test<-morePCs(rld,intgroup=c("condition"),returnData=TRUE,firstPC = 3,secondPC = 4)
# percentVar<-round(100*attr(test,"percentVar"))
# #need to update the axis labels yourself
# ggplot(test, aes(PC1, PC2, color=condition, label=name)) +
#   geom_point(size=3) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   geom_text_repel(aes(label = name),color="black")+
#   xlab(paste0("PC3: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC4: ",percentVar[2],"% variance")) + 
#   coord_fixed()
#get results, adjusted by their correction method
resultsNames(dds)
res<-lfcShrink(dds, coef="condition_disease_vs_control",type="apeglm")

#don't need these variables anymore
rm(pcaData,percentVar,rld,test)

###################################### differentially expressed genes #######################################################
library(tidyverse)
library(data.table)

#add gene names 
gencode_genes<-read.csv("C:/Users/angelawei/Box/_Arboleda Lab Stuff/personal/Angela_Wei/Reference_genome/GencodeV31_GR38/gencode.v31.primary_assembly.genes.csv")
#gencode_genes<-read.csv("/Users/angelawei/Box Sync/_Arboleda Lab Stuff/personal/Angela_Wei/Reference_genome/GencodeV31_GR38/gencode.v31.primary_assembly.genes.csv")
gencode_genes<-dplyr::select(gencode_genes,gene_name,gene_type,seqnames,hgnc_id,gene_id)
hcng_ref<-read.delim("C:/Users/angelawei/Box/_Arboleda Lab Stuff/personal/Angela_Wei/Reference_genome/HCNG/HCNG_05-09-19.txt")
#hcng_ref<-read.delim("/Users/angelawei/Box Sync/_Arboleda Lab Stuff/personal/Angela_Wei/Reference_genome/HCNG/HCNG_05-09-19.txt")
hcng_ref<-dplyr::select(hcng_ref,HGNC.ID,Approved.name,Chromosome,RefSeq.IDs,OMIM.ID.supplied.by.OMIM.)
#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/05_DGE_DESeq2_analysis2.md
#convert results from DESeq2 object to data frame and convert rownames to gene column
res<- res %>%
  data.frame() %>%
  rownames_to_column(var="gene_id") %>%
  as_tibble()

#this dataframe will be one the hold all the extra info (norm counts, HGNC name, chr, etc)
#now add hgnc_id, gene_name, chr,gene_type
res<-merge(res,gencode_genes,all.x=TRUE)
#now add hgnc stuff
res<-merge(res,hcng_ref,all.x=TRUE,by.x="hgnc_id",by.y="HGNC.ID")
#now order by padj
final_res<-res%>%
  filter(padj<0.05)%>%
  #filter(abs(log2FoldChange)>0.6)%>%
  arrange(desc(abs(log2FoldChange)))
#write.csv(final_res,"sigDEGs FC1.5 Issy ASXL1 blood all samples adjGender.csv")
hist(final_res$log2FoldChange)

###################################### gene count plots for DEGs #######################################################
library(ggplot2)
library(ggrepel)
n<-100
for (i in 1:n) {
  geneCounts<-plotCounts(dds,gene = as.character(final_res[i,"gene_id"]),intgroup = "condition",returnData = TRUE)
  geneCounts<-geneCounts%>%
    rownames_to_column(var="sampleName")
  if(is.na(final_res[i,"Approved.name"])==TRUE){
    ggplot(geneCounts, aes(x=condition,y=count,color=condition))+
      geom_point(position=position_jitter(w = 0.1,h = 0),size=3) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle( paste0(final_res[i,"gene_id"], " normalized counts")) +
      geom_text_repel(aes(label = sampleName),color="black")+
      ylab("Normalized counts")+
      labs(subtitle = paste0("log2FC= ",final_res[i,"log2FoldChange"], " , padj=",final_res[i,"padj"]))
    ggsave(file.path(paste0(getwd(),"/normalizedCounts_DEGs/",i,"-",final_res[i,"gene_id"],".png")),
           plot=last_plot(),device = "png",dpi=300)
  }
  else{
    ggplot(geneCounts, aes(x=condition,y=count,color=condition))+
      geom_point(position=position_jitter(w = 0.1,h = 0),size=3) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle( paste0(final_res[i,"gene_name"], ": ",final_res[i,"Approved.name"], " normalized counts")) +
      geom_text_repel(aes(label = sampleName),color="black")+
      ylab("Normalized counts")+
      labs(subtitle = paste0("log2FC= ",final_res[i,"log2FoldChange"], ", padj=",final_res[i,"padj"]))
    ggsave(file.path(paste0(getwd(),"/normalizedCounts_DEGs/",i,"-",gsub("/","_",as.character(final_res[i,"gene_name"])),"-",gsub("/","_",as.character(final_res[i,"Approved.name"])),".png")),
           plot=last_plot(),device = "png",dpi=300)
  }
}

#this function is for extra plots required on a case by case basis
geneCountPlot<-function(geneEnsembl, wantLabels){
  number<-which(res$gene_id==geneEnsembl)
  geneCounts<-plotCounts(dds,gene = geneEnsembl,intgroup = "condition",returnData = TRUE)
  geneCounts<-geneCounts%>%
    rownames_to_column(var="sampleName")
  if(wantLabels==TRUE){
    ggplot(geneCounts, aes(x=condition,y=count))+
      geom_point(position=position_jitter(w = 0.1,h = 0),size=3) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle(paste0(res[number,"gene_name"], " normalized counts"))+
      geom_text_repel(aes(label = sampleName),color="black")+
      ylab("Normalized counts")+
      labs(subtitle = paste0("log2FC= ",round(res[number,"log2FoldChange"],2), " , padj=",formatC(res[number,"padj"], format = "e", digits = 2)))
  }
  else{
    ggplot(geneCounts, aes(x=condition,y=count))+
      geom_point(position=position_jitter(w = 0.1,h = 0),size=3) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggtitle(paste0(res[number,"gene_name"], " normalized counts"))+
      #geom_text_repel(aes(label = sampleName),color="black")+
      ylab("Normalized counts")+
      labs(subtitle = paste0("log2FC= ",round(res[number,"log2FoldChange"],2), " , padj=",formatC(res[number,"padj"], format = "e", digits = 2)))
  }
}

#asxl1
geneCountPlot("ENSG00000171456.19")
#hoxa5 count plot labelled asxl blood rnaseq adjGender
#hoxa5
geneCountPlot("ENSG00000106004.5",wantLabels = FALSE)
#hoxa9
geneCountPlot("ENSG00000078399.18",wantLabels = FALSE)
#hoxb3
geneCountPlot("ENSG00000120093.11",wantLabels = FALSE)
#hoxb4
geneCountPlot("ENSG00000182742.6",wantLabels = FALSE)
#asxl2
geneCountPlot("ENSG00000143970.16",wantLabels = FALSE)
#asxl3
geneCountPlot("ENSG00000141431.12",wantLabels = FALSE)

###################################### heat map for DEGs #######################################################
#now we need the normalized counts
normalized_counts <- counts(dds,normalized=TRUE) %>%
  data.frame() %>%
  rownames_to_column(var="gene_id")
#get the n most sig DEGs by padj
#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
n<-nrow(final_res)
topNDEGs<-final_res %>%
  #arrange(padj) %>%
  arrange(desc(abs(log2FoldChange)))%>%
  pull(gene_id) %>%
  head(n)
#get the normalized counts for the  sig DEGs
topNDEGs_counts<-normalized_counts %>%
  filter(gene_id %in% topNDEGs)%>%
  column_to_rownames(var="gene_id")
meta <- coldata%>%
  rownames_to_column(var="sampleName") %>%
  as_tibble()
#make a heat map based on normalized counts of the top 300 DEGs based on lfc
library(pheatmap)
annotation <- meta %>%
  dplyr::select(sampleName,condition) %>%
  data.frame(row.names = "sampleName")
rownames(annotation)<-colnames(normalized_counts)[2:ncol(normalized_counts)]
#color palette
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlGnBu")
#THIS IS A HEAT MAP OF Z-SCORES OF NORMALIZED COUNTS!!
pheatmap(topNDEGs_counts,
         color = heat_colors,
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#don't need these variables anymore
rm(normalized_counts,n,topNDEGs,topNDEGs_counts,meta,annotation,heat_colors)

###################################### gene ontology #######################################################
#for over-representation analysis
#all genes tested in DE are used for background
#these need to be ensembl IDs
#get rid of version numbers
all_genes<-gsub("\\..*","",gencode_genes$gene_id)
sig_genes<-gsub("\\..*","",final_res$gene_id)
library(clusterProfiler)
#this is the Genome wide annotation for Human
library(org.Hs.eg.db)
ego_all<-enrichGO(gene=sig_genes,
                  universe=all_genes,
                  keyType = "ENSEMBL",
                  OrgDb = org.Hs.eg.db, 
                  ont = "ALL",
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
cluster_summary_all<-data.frame(ego_all)
#write.csv(cluster_summary_all,"geneOntology all DEGs issy asxl1 blood allSamples adjGender.csv")
barplot(ego_all,showCategory = 20)
#simplify similar concept terms
ego_bp<-enrichGO(gene=sig_genes,
                 universe=all_genes,
                 keyType = "ENSEMBL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.05, 
                 readable = TRUE)
#https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
ego_bp_simp<-clusterProfiler::simplify(x=ego_bp,cutoff=0.7,by="p.adjust",select_fun=min)
summary_simp<-data.frame(ego_bp_simp)
emapplot(ego_bp_simp)
barplot(ego_bp_simp,showCategory = 20)
#write.csv(summary,"simplified geneOntology all DEGs issy asxl1 blood allSamples adjGender.csv")

# barplot(ego_bp,showCategory = 50)
# dotplot(ego_bp,showCategory=20)
# emapplot(ego_bp)
# rm(all_genes,sig_genes,ego_all,cluster_summary_all,ego_bp)

########for homer motifs######
final_res$ensembl_noVersion<-gsub("\\..*","",final_res$gene_id)
rnaTxt<-data.frame(final_res[,15])
#write.table(rnaTxt,"all_DEGs_issy_asxl1_blood.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)
up<-final_res[final_res$log2FoldChange>0,]
rnaTxt<-data.frame(up[,15])
#write.table(rnaTxt,"up_DEGs_issy_asxl1_blood.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)
down<-final_res[final_res$log2FoldChange<0,]
rnaTxt<-data.frame(down[,15])
#write.table(rnaTxt,"down_DEGs_issy_asxl1_blood.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)

#########################volcano plot###############################################
library(tidyverse)
#column for if peak meets lfc and/or padj threshold
res<-res %>% mutate(threshold = padj<0.05 & abs(log2FoldChange)>=0.58)
#order everything from increasing padj and add a blank column to hold approved name
res<-res %>% arrange(padj) %>% mutate(genelabels = "")
#fill the top 30 genes with an approved name
res$genelabels[1:30]<-as.character(res$gene_name[1:30])
title<-"Differential expression, ASXL1 BO Blood vs controls"
#this is a version of the volcano plot without gene names
ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=factor(threshold)))+
  geom_point()+
  ggtitle(title) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme_bw()+
  geom_hline(yintercept=1.301, linetype="dashed", 
             color = "grey", size=1)+
  geom_vline(xintercept=0.58, linetype="dashed", 
             color = "grey", size=1)+
  geom_vline(xintercept=-0.58, linetype="dashed", 
             color = "grey", size=1)+
  scale_color_manual(breaks = c("FALSE", "TRUE"),
                     values=c("#0066CC", "#FF0000"))+
  labs(color="padj<0.05 & abs(log2FC)>0.58")
#if you get error: "Warning message: Removed 34954 rows containing missing values (geom_point)." it is because outliers will not have a padj
#this is a version of the volcano plot where top 30 DEGs (by padj) are labeled with their names
library(ggrepel)
ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=factor(threshold)))+
  geom_point()+
  ggtitle(title) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme_bw()+
  geom_hline(yintercept=1.301, linetype="dashed", 
             color = "grey", size=1)+
  geom_vline(xintercept=0.58, linetype="dashed", 
             color = "grey", size=1)+
  geom_vline(xintercept=-0.58, linetype="dashed", 
             color = "grey", size=1)+
  scale_color_manual(breaks = c("FALSE", "TRUE"),
                     values=c("#0066CC", "#FF0000"))+
  labs(color="padj<0.05 & abs(log2FC)>0.58")+
  geom_text_repel(aes(label = genelabels),color="black")
