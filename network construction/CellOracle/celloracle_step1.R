##############################################
#####             CellOracle             #####
##### 1. scATAC-seq analysis with Cicero #####
##############################################
library(cicero)
library(data.table)
library(monocle3)
library(tidyr)

# Load data
data_folder <- "/home/fuyq/GRN/single_cell/data/pbmc/matrix"
indata <- Matrix::readMM(paste0(data_folder, "/matrix.mtx"))
indata@x[indata@x > 0] <- 1

cellinfo <- read.table(paste0(data_folder, "/barcodes.tsv"))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

featureinfo <- fread(paste0(data_folder, "/features.tsv"))
featureinfo<-as.data.frame(featureinfo)

row.names(indata) <- featureinfo$V2
colnames(indata) <- row.names(cellinfo)

rows_with_colon <- grepl(":", rownames(indata))
table(rows_with_colon)
indata <- indata[rows_with_colon, ]
rownames(indata) <- gsub("[:-]", "_", rownames(indata))

peak <- data.frame(site_name=rownames(indata))
peak <- peak %>%
  separate(site_name, into = c("chr", "bp1", "bp2"), sep = "_", convert = TRUE,remove = FALSE)
peakinfo<-peak[,c(2,3,4,1)]
row.names(peakinfo) <- peakinfo$site_name

# Create monocle3 CDS object
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

hist(Matrix::colSums(exprs(input_cds)))
max_count <-  20000
min_count <- 1000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count]
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count]

# Dimensionality reduction (LSI + UMAP)
set.seed(2026)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# Run Cicero (co-accessibility analysis)
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
save_monocle_objects(cds = cicero_cds,
                     directory = "/home/fuyq/GRN/single_cell/celloracle/pbmc10k/result/cellorcale/BaseGRN_input_data_preparation/monocle_save")
chromosome_length <- read.table("/home/fuyq/GRN/single_cell/celloracle/chr_length/hg38.chrom.sizes.txt",header = F)
conns <- run_cicero(cicero_cds, chromosome_length)
saveRDS(conns,file ="/home/fuyq/GRN/single_cell/celloracle/pbmc10k/result/cellorcale/BaseGRN_input_data_preparation/run_cicero.RDS")
# Export peaks and Cicero results
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = "/home/fuyq/GRN/single_cell/celloracle/pbmc10k/result/cellorcale/BaseGRN_input_data_preparation/all_peaks.csv")
write.csv(x = conns, file = "/home/fuyq/GRN/single_cell/celloracle/pbmc10k/result/cellorcale/BaseGRN_input_data_preparation/cicero_connections.csv")

