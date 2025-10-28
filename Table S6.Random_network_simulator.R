library(igraph)

# Set the size of random network
n=162   #number of nodes
e=1132  #number of edges

n=31   #number of nodes
e=72  #number of edges

n=32   #number of nodes
e=124  #number of edges

n=60   #number of nodes
e=333  #number of edges
# Read otu.txt (make sure this file exists in the working directory)

# Generate 1000 random networks
for (i in 1:1000) {
  
  # Generate a random network
  g <- erdos.renyi.game(n, e, 'gnm', directed = FALSE)
  
  # Convert the graph to an adjacency matrix
  adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

  
  # Global topological features
  c <- cluster_walktrap(g)
  md <- modularity(g, membership(c), weights = NULL)
  cc <- transitivity(g, vids = NULL, weights = NULL)
  apl <- average.path.length(g, directed = FALSE, unconnected = TRUE)
  gd  <- graph.density(g, loops = FALSE)
  nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
  
  node.degree <- degree(g, v = V(g), mode = "all")
  ad  <- mean(node.degree)
  
  global.topol <- data.frame(n, e, cc, apl, md, gd, nd, ad)
  
  # Save global topological features to file
  write.table(global.topol, file = sprintf("random_network.xls", n, e),
              
              append = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
