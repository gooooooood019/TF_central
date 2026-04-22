##############################################
#####              构建网络              #####
#####           scTenifoldKnk            #####
##############################################
library(dplyr)
library(scTenifoldKnk)
library(ggplot2)
library(Matrix)
library(data.table)
library(parallel)
# Load input data
WT1 <- readMM('/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/data/example/matrix_WT.mtx')
rownames(WT1) <- readLines('/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/data/example/geneNames_WT.txt')
colnames(WT1) <- readLines('/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/data/example/barcodes_WT.txt')

# Quality control and filtering
source('/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/code/scQC.R')
WT1 <- scQC(WT1)
WT1 <- WT1[rowMeans(WT1 != 0) > 0.05,]
WT1 <- WT1[!grepl('^Rp[[:digit:]]+|^Rpl|^Rps|^Mt-', rownames(WT1), ignore.case = TRUE),]

# Define candidate genes
tf<-fread("/home/fuyq/GRN/GTEx/tf_mm_new.txt",header=F)
Gene<-intersect(rownames(WT1),tf$V1)

# Virtual knockout analysis (scTenifoldKnk)
outdir <- "/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/result/example/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (gene in Gene) {
  set.seed(1)
  cat("Processing:", gene, "\n")
  
  WT <- WT1
  result <- tryCatch({
    scTenifoldKnk(WT, gKO = gene)
  }, error = function(e) {
    cat("Error occurred for gene", gene, ":", e$message, "\n")
    return(NULL)
  })
  
  # Save results if successful
  if (!is.null(result)) {
    file_name <- paste0(gene, "_GSE130626.rds")
    saveRDS(result, file = file.path(outdir, file_name))
    cat("Saved results for", gene, "as", file_name, "\n")
  }
}

# Extract significant differential regulation
files<-list.files("/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/result/example",full.names = T)
files_names<-list.files("/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/result/example")
files_names<-sub("_GSE130626.rds","",files_names)
numCores<-32
result_list <- mclapply(seq_along(files), function(i) {
  tmp <- readRDS(files[i])
  tmp_diff<-tmp$diffRegulation
  tmp_diff<-tmp_diff[tmp_diff$p.value<0.05,]#p.adj<0.05
}, mc.cores = numCores)
names(result_list) <- files_names
saveRDS(result_list,file="/home/fuyq/GRN_process/virtual_knock/scTenifoldKnk/result/result_list_v2.RDS")
