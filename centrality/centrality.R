##############################################
#####             计算中心性             #####
##############################################
library(data.table)
library(igraph)
library(centiserve)

# Load adjacency list (GRN results from SCENIC)
adj_ctx<-readRDS("/home/fuyq/GRN/single_cell/scenic/PBMC10k/centrality/adj_ctx_list.RDS")
# Full list of structural factors (SFs)
SF<-c("CTCF", "RAD21","SMC1A","SMC3","STAG1","STAG2","NIPBL","MAU2","WAPL","PDS5A","PDS5B",
      "ESCO1","ESCO2","HDAC8","CDCA5","YY1","ZNF143","HDGF","CGGBP1","ADNP","HMGB2","GATAD2A",
      "ZNF512","ZNF532","MORC2","GATAD2B","MORC3","GABPA","NRF1","ZBTB11","JARID2","ZNF281",
      "MGA","PRDM10","ZNF296","NR0B1","POU5F1","ZNF462","ESRRB","DPPA2") #40
# the list of common structural factors (cSFs)
cSF<-c("CTCF", "RAD21","SMC1A","SMC3","STAG1","STAG2","NIPBL","MAU2","WAPL","PDS5A","PDS5B",
       "ESCO1","ESCO2","HDAC8","CDCA5","YY1","ZNF143") #17
# the list of uncommon structural factors (uSFs)
uSF<-c("HDGF","CGGBP1","ADNP","HMGB2","GATAD2A","ZNF512","ZNF532","MORC2","GATAD2B","MORC3","GABPA","NRF1",
       "ZBTB11","JARID2","ZNF281","MGA","PRDM10","ZNF296","NR0B1","POU5F1","ZNF462","ESRRB","DPPA2") #23
#####************direct graph************#####
results <- data.frame()
for (i in seq_along(adj_ctx)) {
  
  ### Load data
  df<-adj_ctx[[i]]
  # Identify all regulators in the network
  tf<-unique(df$TF)
  # Classify nodes into SFs and TFs
  SF_new<-intersect(tf,SF)
  TF_new <-setdiff(tf,SF_new)
  
  ### Construct directed graph
  graph<-graph_from_data_frame(df, directed = TRUE)
  graph <- igraph::simplify(graph, remove.multiple = F, remove.loops = TRUE)
  
  ### Centrality calculations
  # 1. degree centrality
  SF_degree<-igraph::degree(graph,mode = "out",v = SF_new)
  TF_degree<-igraph::degree(graph,mode = "out",v = TF_new)
  # 2. Betweenness centrality
  SF_betweenness<-betweenness(graph,v = SF_new)
  TF_betweenness<-betweenness(graph,v = TF_new)
  # 3. Closeness centrality
  SF_closeness<-closeness(graph, mode = "out",v = SF_new)
  TF_closeness<-closeness(graph, mode = "out",v = TF_new)
  # 4. Eigenvector centrality
  eigenvector<-eigen_centrality(graph)$vector
  SF_eigen<-eigenvector[SF_new]
  TF_eigen<-eigenvector[TF_new]
  # 5. PageRank centrality
  pagerank_scores<-page_rank(graph)$vector
  pagerank<-data.frame(Gene = names(pagerank_scores), PageRank = pagerank_scores)
  pagerank$Category <- ifelse(pagerank$Gene %in% SF_new, "SF",
                              ifelse(pagerank$Gene %in% TF_new, "NSF", "Other"))
  SF_pagerank <- pagerank$PageRank[pagerank$Category == "SF"]
  TF_pagerank <- pagerank$PageRank[pagerank$Category == "NSF"]
  # 6. Radiality centrality
  rad<-radiality(graph)
  SF_rad <- rad[SF_new]
  TF_rad <- rad[TF_new]
  # 7. Local clustering coefficient
  local_cc <- transitivity(graph, type = "local")
  SF_cc <- local_cc[SF_new]
  TF_cc <- local_cc[TF_new]
  
  ### Store results
  results<-rbind(results, data.frame(
    degree = c(SF_degree,TF_degree),
    betweenness = c(SF_betweenness,TF_betweenness),
    closeness = c(SF_closeness,TF_closeness),
    eigenvector = c(SF_eigen,TF_eigen),
    pagerank = c(SF_pagerank,TF_pagerank),
    radiality = c(SF_rad, TF_rad),
    Clustering_coefficient = c(SF_cc, TF_cc),
    tissue = rep(names(adj_ctx)[i],length(unique(df$TF))),
    type = c(rep("SF", length(SF_degree)),rep("TF", length(TF_degree))),
    Gene = c(names(SF_degree),names(TF_degree))
  ))
}
  ### cSF vs uSF vs TF
results$type3 <- ifelse(results$Gene %in% cSF, "cSF",
                        ifelse(results$Gene %in% uSF, "uSF", "TF"))

saveRDS(results,file = "/home/fuyq/GRN/single_cell/scenic/PBMC10k/centrality/direct/result_ctx.rds")
#####************undirect graph************#####
results_list <- list()
for (i in seq_along(adj_ctx)) {
  
  tissue_name <- names(adj_ctx)[i]
  cat("Processing tissue:", tissue_name, "\n")
  
  ### Load data
  df <- adj_ctx[[i]]
  # Identify all regulators in the network
  tf <- unique(df$TF)
  # Classify nodes into SFs, TFs and Genes
  SF_new <- intersect(tf, SF)
  TF_new <- setdiff(tf, SF_new)
  gene <- unique(df$Target)
  Gene_new <- setdiff(gene, tf)
  
  ### Construct undirected graph
  graph <- graph_from_data_frame(df, directed = FALSE)
  graph <- igraph::simplify(graph, remove.multiple = F, remove.loops = TRUE)
  
  ### Centrality calculations
  # 1. degree centrality
  degree <- degree(graph)
  SF_degree <- degree[SF_new]
  TF_degree <- degree[TF_new]
  Gene_degree <- degree[Gene_new]
  # 2. Betweenness centrality
  betweenness <- betweenness(graph)
  SF_betweenness <- betweenness[SF_new]
  TF_betweenness <- betweenness[TF_new]
  Gene_betweenness <- betweenness[Gene_new]
  # 3. Closeness centrality
  closeness <- closeness(graph)
  SF_closeness <- closeness[SF_new]
  TF_closeness <- closeness[TF_new]
  Gene_closeness <- closeness[Gene_new]
  # 4. Eigenvector centrality
  eigenvector <- eigen_centrality(graph)$vector
  SF_eigen <- eigenvector[SF_new]
  TF_eigen <- eigenvector[TF_new]
  Gene_eigen <- eigenvector[Gene_new]
  # 5. PageRank centrality
  pagerank_scores <- page_rank(graph)$vector
  pagerank <- data.frame(Gene = names(pagerank_scores), PageRank = pagerank_scores)
  pagerank$Category <- ifelse(pagerank$Gene %in% SF_new, "SF",
                              ifelse(pagerank$Gene %in% TF_new, "TF", "Gene"))
  
  SF_pagerank <- pagerank$PageRank[pagerank$Category == "SF"]
  TF_pagerank <- pagerank$PageRank[pagerank$Category == "TF"]
  Gene_pagerank <- pagerank$PageRank[pagerank$Category == "Gene"]
  # 6. Radiality centrality
  rad<-radiality(graph)
  SF_rad <- rad[SF_new]
  TF_rad <- rad[TF_new]
  Gene_rad <- rad[Gene_new]
  # 7. Local clustering coefficient
  local_cc <- transitivity(graph, type = "local")
  SF_cc <- local_cc[SF_new]
  TF_cc <- local_cc[TF_new]
  Gene_cc <- local_cc[Gene_new]
  
  ### Store results
  results <- data.frame(
    degree = c(SF_degree, TF_degree, Gene_degree),
    betweenness = c(SF_betweenness, TF_betweenness, Gene_betweenness),
    closeness = c(SF_closeness, TF_closeness, Gene_closeness),
    eigenvector = c(SF_eigen, TF_eigen, Gene_eigen),
    pagerank = c(SF_pagerank, TF_pagerank, Gene_pagerank),
    radiality = c(SF_rad, TF_rad, Gene_rad),
    Clustering_coefficient = c(SF_cc, TF_cc, Gene_cc),
    tissue = rep(names(adj_ctx)[i], length(c(SF_degree, TF_degree, Gene_degree))),
    type = c(rep("SF", length(SF_degree)),
             rep("TF", length(TF_degree)),
             rep("Gene", length(Gene_degree))),
    Gene = c(names(SF_degree),names(TF_degree),names(Gene_degree))         
  )
  results_list[[names(adj_ctx)[i]]] <- results
}

results_list <- lapply(results_list, function(df){
  df$type3 <- df$type
  df$type3[df$Gene %in% uSF] <- "uSF"
  df$type3[df$Gene %in% cSF] <- "cSF"
  return(df)
})
saveRDS(results_list,file = "/home/fuyq/GRN/single_cell/scenic/PBMC10k/centrality/undirect/result_ctx.rds")

