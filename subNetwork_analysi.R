
################## OTU filtering, network generation, topological analysis and export OTU table ###############################
gill_gn1<-g



T_list <- list.files(path="F:/net/2",
                     pattern="*otu.txt")

V.gill <- list()
for (n in T_list){
  n.fish <-  read.table(n)
  OTU_table <- n.fish;
  table <- OTU_table
  table[table > 1]<-1;
  
  #Filter out empty rows
  dt<- OTU_table
  dt[dt>1]<-1
  no<-which(rowSums(dt)/ncol(dt)>0.2)
  OTU_Abu <- OTU_table[no,]
  #node name
  V.gill[[n]] <-row.names(OTU_Abu)
}

for (i in V.gill){  
  V.sub = subset(i, i %in% V(gill_gn1)$name) 
  g <- induced.subgraph(gill_gn1, vids = V.sub, impl = c("copy_and_delete")) 
  g1 <-simplify(g)
  g2 <-delete.vertices(g1, names(degree(g1)[degree(g1)==0]))
  V(g2)$label <- V(g2)$name
  
  ### Calculating subnetwork topological properties
  c <- cluster_walktrap(g2)
  # Global toplogical features
  modularity(c)
  md <- modularity(g2, membership(c), weights = NULL)
  cc <- transitivity(g2, vids = NULL,
                     weights = NULL)
  apl <- average.path.length(g2, directed=FALSE, unconnected=TRUE)
  gd  <- graph.density(g2, loops=FALSE)
  nd  <- diameter(g2, directed = FALSE, unconnected = TRUE, weights = NULL)
  node.degree <- degree(g2, v = V(g2), mode="all")
  ad  <- mean(node.degree)
  e <- ecount(g2)
  v <- vcount(g2)
  global.topology <- data.frame(e,v,cc,apl,md,gd,nd,ad)
  write.csv(global.topology, paste0("global.topology", length(i), ".csv"))
  
  # Node toplogical features
  betweenness.centrality = centralization.betweenness(g2)$centralization 
  
  closeness.centrality <- closeness(g2, vids = V(g2),
                                    weights = NULL, normalized = FALSE)
  node.transitivity <- transitivity(g2, type = c("local"), vids = NULL,
                                    weights = NULL)
  
  node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
  write.csv(node.topology, paste0("node.topology", length(i), ".csv"))
  write.graph(g2, paste0("network", length(i), ".gml"), format='gml')}

