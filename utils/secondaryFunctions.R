# function to create links
createLink <- function(val) {
  sprintf('<a href="https://reactome.org/PathwayBrowser/#/%s" target="_blank" class="btn btn-primary">more info</a>',val)
}

# network descriptives
netsummary <- function(net) {
  nodes = vcount(net)
  edges = ecount(net)
  diameter = max(eccentricity(net, mode = "all"))
  comps = components(net)$no
  shortestpath = round(mean_distance(net, directed = F), 3)
  dens = round(graph.density(net), 3)
  av_neigh = round(mean(degree(net, v=V(net), mode = "all", loops = F)), 3)
  cc = round(transitivity(net, type = "average"), 3)
  centr = round(centr_degree(net)$centralization, 3)
  return(data.frame(Nodes = nodes, 
                    Edges = edges, 
                    Diameter = diameter, 
                    Shortest_Path = shortestpath,
                    Density = dens,
                    Average_Neighbors = av_neigh,
                    Clustering_Coefficient = cc,
                    Centralization = centr,
                    Components = comps))
}

# network topologies
toposummary <- function(net) {
  cc <- transitivity(net, type = "local")
  cc[is.na(cc)] <- 0
  deg <- degree(net, v=V(net), mode = "all", loops = F, normalized = F)
  deg_norm <- degree(net, v=V(net), mode = "all", loops = F, normalized = T)
  clo <- closeness(net, vids = V(net), mode = "all", normalized = F)
  bet <- betweenness(net, v=V(net), directed = F, normalized = T)
  eig <- eigen_centrality(net, directed = F, scale = F)
  pr <- page_rank(net, vids = V(net), directed = F, weights = E(net)$weight)
  ecc <- eccentricity(net, vids = V(net), mode = c("all"))
  edgebet <- edge_betweenness(net, e = E(net), directed = FALSE)
  weights <- E(net)$weight
  node_at <- data.frame(degree = deg,
                        norm_degree = deg_norm,
                        clust_coef = cc,
                        cc = round(cc, 2),
                        closeness = clo,
                        betwenness = bet,
                        eigencent = eig$vector,
                        eccentricity = ecc,
                        pagerank = pr$vector,
                        pr = round(pr$vector, 3),
                        labels = names(V(net)))
  return(node_at)
}


remDownComp <- function(edgeList, nodeList) {

  comps <- components(graph_from_data_frame(edgeList, directed = T))$membership
  comps <- data.frame(cluster = comps, class = nodeList$fun[match(names(comps), nodeList$id)])
  comps2rem <- table(comps$cluster, comps$class)
  comps2rem <- as.integer(row.names(comps2rem)[comps2rem[,"Receptor"]!=0])
  comps <- comps[comps$cluster%in%comps2rem, ]
  ord <- as.numeric(as.character(names(sort(table(comps$cluster), decreasing = T))))
  pos <- 1:length(ord)
  comps$cluster <- pos[match(unname(comps$cluster), ord)]
  edgeList <- edgeList[edgeList$from%in%row.names(comps) | edgeList$to%in%row.names(comps), ]
  nodeList$component <- comps$cluster[match(nodeList$id, row.names(comps))]
  nodeList <- nodeList[!is.na(nodeList$component) & nodeList$id%in%unique(c(edgeList$from, edgeList$to)), ]
  return(list(edgeList, nodeList))
}

sortComps <- function(nodeList) {
  ord <- as.numeric(as.character(names(sort(table(nodeList$components), decreasing = T))))
  pos <- 1:length(ord)
  nodeList$components <- pos[match(unname(nodeList$components), ord)]
  return(nodeList)
}

changeColorOfNodes <- function(nodes, selected.nodes) {
  #tmp <- mask[match(selected.nodes, mask$id), ]
  #nodes[match(tmp$id, nodes$id), c(1,9)] <- tmp[, c(1,3)]
  nodes %>% mutate(color = if_else(id %in% selected.nodes, "gray25", color))
}

warnColumn <- HTML(paste('<br/>', '<br/>', "Check your input data, lacks list specification column"))
warnEntries <- HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists"))
warnMaxlist <- HTML(paste('<br/>', '<br/>', "Check your input data, No more than 8 different lists supported!"))
warnQuenya <- HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!"))
warnCommon <- HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!"))
warnFormat <- HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed"))
warnAnnot <- HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!"))
warnInt <- HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!"))
warnScore <- HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -"))
warnCrosst <- HTML(paste('<br/>', '<br/>', "No interactions found!"))
warnDrug <- HTML(paste('<br/>', '<br/>', "No Drug interactions found with Receptors provided"))
warnEnr <- HTML(paste('<br/>', '<br/>', "No terms selected for enrichment!"))

  
