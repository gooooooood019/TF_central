##############################################
#####              构建网络               #####
#####             DIRECT-NET             #####
##############################################
library(DIRECTNET)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2016)
library(presto)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1000 * 1024^2)

pbmc <- readRDS("/home/fuyq/GRN/single_cell/saveRDS/pbmc_wnn_analysis.Rds")

# Perform differential expression using Wilcoxon AUC test
markers_all <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_all <- markers_all[which(markers_all$auc > 0.5), , drop = FALSE]
markers <- data.frame(gene = markers_all$feature,
				       group = markers_all$group)\
# Create marker list per cell type
celltypes <- unique(markers$group)
marker_list <- list()
for (i in 1:length(celltypes)) {
  marker1 <- markers_all[markers$group == celltypes[i], ]
  marker_list[[celltypes[i]]] <- as.character(marker1$feature[marker1$auc > 0.5])
}

markers_groups <- unique(unlist(marker_list)) #
markers_groups <- lapply(markers_groups, function(x) strsplit(x,"[.]")[[1]][1])
markers_groups <- unique(unlist(markers_groups))

# Run DIRECT-NET on ATAC object
pbmc_atac <- Run_DIRECT_NET(pbmc_atac, peakcalling = FALSE,
                       k_neigh = 50, atacbinary = TRUE,
                       max_overlap=0.5, size_factor_normalize = FALSE,
                       genome.info = genome_info, focus_markers = markers_groups )
saveRDS(pbmc_atac,file = "/home/fuyq/GRN/single_cell/direct_net/PBMC10k/direct_net/PBMC_direct.net.RDS")

direct.net_result <- Misc(pbmc_atac, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind,direct.net_result)) # links for markers
direct.net_result$function_type <- gsub("HF","HC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("Rest","MC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("LF","LC",direct.net_result$function_type)

# Identify Differential Accessible Peaks (DA peaks)
groups <- unique(pbmc@meta.data$celltype)

da_peaks_list <- list()
for (i in seq_along(groups)) {
  message("Processing cell type: ", groups[i])
 
  da_peaks <- FindMarkers(
    object = pbmc,
    ident.1 = groups[i],
    group.by = "celltype",
    min.pct = 0.2,
    logfc.threshold = 0.6,
    test.use = 'LR',
    only.pos = TRUE
  )
  da_peaks_list[[groups[i]]] <- da_peaks
}

# Generate CRE to gene links
CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = markers)

# Identify CREs overlapping differential peaks
Focused_CREs <- generate_CRE(L_G_record = CREs_Gene$distal, P_L_G_record = CREs_Gene$promoter, da_peaks_list)

# Detect TF binding on distal CREs
L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$distal, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = markers)

# Detect TF binding on promoter regions
P_L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$promoter, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = markers)

# Generate TF-gene connections
network_links <- generate_links_for_Cytoscape(L_G_record = Focused_CREs$L_G_record, L_TF_record, P_L_G_record = Focused_CREs$P_L_G_record, P_L_TF_record,groups)
saveRDS(network_links, file = "/home/fuyq/GRN/single_cell/direct_net/PBMC10k/direct_net/network_links.RDS")
