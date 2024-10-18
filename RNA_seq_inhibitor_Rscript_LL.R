#R script of code, (I was having trouble with R markdown)
#setting working directory to lab box folder
setwd("C:\\Users\\leall\\Box\\_Arboleda Lab Stuff\\personal\\Lesly_Leal\\RNA_seq_KAT6Ainhibitor_small_dataset")

library(tidyverse)
library(readxl)
library(dplyr)
library(pkgbuild) ##load pkgbuild
library(DESeq2) ##Load in Deseq library


countdata <-read_excel("C:/Users/leall/OneDrive/KAT6A Inhibitor small dataset/KAT6A-inhibitors-controlNPC-CountMatrix_AN008D_221228 (1).xlsx")
coldata <- data.frame(read_excel("C://Users//leall//Box//_Arboleda Lab Stuff//personal//Lesly_Leal//KAT6A Inhibitor small dataset//KAT6A-inhibitors-controlNPC-colData-sampleINFO_AN008D_221228.xlsx"))


# Step 1: Convert tibble to data frame (if it's a tibble)
countdata_df <- as.data.frame(countdata)

# Step 2: Set row names to the first column (gene IDs)
#rownames(countdata_df) <- countdata_df[, 1]  # Set row names to the first column
#countdata_df <- countdata_df[, -1]            # Remove the first column


filtered_coldata <- coldata %>%
  filter((timept == "24H"))
#view(filtered_coldata)

#extract sample IDs
rownames(coldata) <- coldata[, 2]  # Here, `coldata[, 2]` refers to the second column

sampleIDs <- filtered_coldata$sampleID

# Step 2: Convert the count data to a data frame if it's not already done
countdata_df <- as.data.frame(countdata)

# Step 3: Identify the gene ID column (assuming it's the first column)
geneID_column <- colnames(countdata_df)[1] 

# Step 4: Filter the count data to only include samples in the filtered coldata
filtered_countdata <- countdata_df %>%
  select(geneID_column, all_of(sampleIDs))

# Step 5: Check the filtered count data
print(head(filtered_countdata))

#reset Variables to coldata and countdata

countdata<- filtered_countdata
coldata <- filtered_coldata

view(countdata)
view(coldata)

##Align rows of coldata co columns of countdata
print(rownames(coldata))
print(coldata[,2])
rownames(coldata) <- coldata[,2]
coldata <- coldata[,-c(1,2)]
view(coldata)

print(rownames(countdata))
print(countdata[,1])
rownames(countdata) <- countdata[,1]
countdata<- countdata[,-1]
view(countdata)
###############################################TRYING TO FILTER OUT DATA WITHOUT COLLAPSING MATRIX####################

# Set row names if your first column is gene names (adjust column index if needed)
#rownames(countdata) <- countdata[[1]]

# Convert remaining data to a numeric matrix
countdata <- as.matrix(countdata)

# Replace all negative values with 0, while keeping 0s and positive values unchanged
#countdata[countdata < 0] <- 0

print(dim(countdata))  

print(colnames(countdata))  # This should give you the sample names in countdata
print(rownames(coldata))    # This should give you the sample names in coldata

# Align column names and row names, if necessary
rownames(coldata) <- colnames(countdata)

print(colnames(countdata)) #CHECKING AGAIN if colnames oc countdata match rownames of col data
print(rownames(coldata)) 

all(colnames(countdata)== rownames(coldata))
class(countdata)

countdata <- round(countdata)  # Round any decimal values to nearest integer
mode(countdata) <- "integer"   # Explicitly set mode to integer

####(running Deseq)##

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition )
##more filtering
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
dds


#set factor level
dds$contition <-relevel(dds$condition, ref= "untreated")

#collapse technical replicates (ASK AILEEN IF THIS IS NECESSARY)

#run DEseq 
dds <- DESeq(dds)

res<- results(dds)
res

#CONTRAST??


#MAplot
plotMA(res)
