### August 1, 2023 ###
# load packages
#devtools::install_github("satijalab/seurat", "seurat5")
library(Seurat)
library(SeuratObject)
library(tidyverse)
#library(Matrix)
library(data.table)
#devtools::install_github("satijalab/azimuth", "seurat5")
library(Azimuth)
#library(biomaRt)
#library(org.Hs.eg.db)
#install.packages("genio")
library(genio)


## read in arguments following terminal command
# order of arguments: sampleID.txt input_dir_counts_matrix output_dir output_name.rds
args <- commandArgs(trailingOnly = TRUE)

############################GEX PROCESSING####################################################

### make gex matrix list and rename col names
gex.list <- list()
gex.seuratv5.list <- list()
n_pcr_samples <- count_lines(args[1])
sample_names <- read.table(args[1])

# if converting to ensemble times out you will have to start the for loop over again
# but can change the starting number from 1 to pick up on the last object made
for (i in 1:n_pcr_samples) {
  # read in counts
  name <- sample_names$V1[i]
  c <- fread(file = 
               paste0(args[2],"/",name, "_counts.tsv.gz"))
  # save df without gene column and set gene column as rownames
  rownames(c) <- c$gene
  m <- as.matrix(c[,2:ncol(c)], nrow = nrow(c))
  rownames(m) <- c$gene
  rownames(m) <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", rownames(m))
  n <- Azimuth:::ConvertEnsembleToSymbol(m, species = "human")
  gex.list[[i]] <- n
  # saving seurat objects
  obj <- CreateSeuratObject(counts = n, min.features = 200)
  obj$pcr_barcode <- name
  gex.seuratv5.list[[i]] <- obj
}

# combining all seurat objects together
if (length(gex.seuratv5.list) > 1) {
  objv5 <- merge(gex.seuratv5.list[[1]],
                 gex.seuratv5.list[2:length(gex.seuratv5.list)],
                 add.cell.ids = sample_names$V1)
} else {
  objv5 <- gex.seuratv5.list[[1]]
}

objv5$sample <- gsub("_[0-9]","",objv5$pcr_barcode)


cells <- colnames(objv5)
rtlig.barcodes <- str_split_fixed(cells, "[0-9]_", n=2)[,2]
ligation.barcodes <-  substr( rtlig.barcodes , start = 1, stop = 10)
rt.barcodes <- substr( rtlig.barcodes , start = 11 , stop = 20)

objv5$rt.barcode <- rt.barcodes
objv5$lig.barcode <- ligation.barcodes

saveRDS(objv5,paste0(args[3],"/",args[4]))