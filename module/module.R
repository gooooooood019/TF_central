##############################################
#####             网络模块化             #####
##############################################
library(igraph)

### Load data
adj_ctx<-readRDS("/home/fuyq/GRN/single_cell/scenic/PBMC10k/centrality/adj_ctx_list.RDS")
outdir<-"/home/fuyq/GRN/single_cell/scenic/PBMC10k/module/csv/"

### Community detection
set.seed(2025)
tissue<-names(adj_ctx)
community_list <- list()
for (i in seq_along(adj_ctx)) {
  tmp<-adj_ctx[[i]]
  # Construct undirected graph
  graph<-graph_from_data_frame(tmp, directed = F)
  # Perform community detection using the Louvain algorithm
  community <- cluster_louvain(graph)
  # Extract node-to-community assignments
  community_data <- data.frame(node = V(graph)$name, community_id = community$membership)
  # Save community assignment for current network
  write.csv(community_data, 
            file =paste0(outdir,"community_",tissue[i],".csv"), 
            row.names = FALSE)
  # Store results in list
  community_list[[i]] <- community_data
}
names(community_list)<-tissue
saveRDS(community_list,file = "/home/fuyq/GRN/single_cell/scenic/PBMC10k/module/module_community.RDS")
