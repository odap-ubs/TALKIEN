# function to create links
createLink <- function(val) {
  sprintf('<a href="https://reactome.org/PathwayBrowser/#/%s" target="_blank" class="btn btn-primary">more info</a>',val)
}

# permutation test distribution
expInts <- function(dat, graph, ddbb, annot) {
 
  nCond <- unique(dat)
  lists <- combn(nCond, m = 2, simplify = T)
  lists <- apply(lists, 2, function(x) paste(x, collapse = " X "))
  res <- data.frame(matrix(NA, nrow = 1000, ncol = length(lists)+1))
  names(res) <- c(lists, "total")
  for (i in 1:1000) {
    nodeList <- list()
    for (j in 1:length(nCond)) {
      nodeList[[j]] <- sample(annot$gene.symbol, replace = F, size = table(dat)[[j]])
    }
    names(nodeList) <- nCond
    nodeList <- lapply(nodeList, function(cond) { 
      ind <- cond%in%row.names(ddbb)
      id <- cond[ind]
      class <- ddbb$class[match(id, row.names(ddbb))]
      cond <- cbind.data.frame(id, class)
      return(cond)})
    nodeList <- lapply(seq_along(nodeList), function(i) {
      list <- rep(names(nodeList)[[i]], nrow(nodeList[[i]]))
      return(cbind.data.frame(nodeList[[i]], list))
    })
    nodes <- unlist(lapply(nodeList, function(x) x$id))
    nodes <- names(table(nodes)[table(nodes)==1])
    nodeList <- lapply(nodeList, function(cond) {cond <- cond[cond$id%in%nodes, ]; return(cond)})
    nodes <- do.call(rbind.data.frame, nodeList)
    nodes$stringId <- annot$string.id[match(nodes$id, annot$gene.symbol)]
    nodeIds <- names(V(graph))[names(V(graph))%in%nodes$stringId]
    nodes <- nodes[nodes$stringId%in%nodeIds, ]
    net <- induced_subgraph(as.undirected(graph), nodeIds)
    net <- set.vertex.attribute(net, "class", value = nodes$class[match(names(V(net)), nodes$stringId)]) 
    net <- set.vertex.attribute(net, "list", value = nodes$list[match(names(V(net)), nodes$stringId)])
    edges <- as.data.frame(as_edgelist(net))
    from <- nodes[match(edges[,1], nodes$stringId), c("list", "class")]
    to <- nodes[match(edges[,2], nodes$stringId), c("list", "class")]
    edges <- cbind.data.frame(edges, from, to)
    names(edges) <- c("from", "to", "list_from", "class_from", "list_to", "class_to")
    
    numEdges <- sapply(lists, function(x) paste(edges$list_from, edges$list_to, sep = " X ") %in% x | paste(edges$list_to, edges$list_from, sep = " X ") %in% x)
    if (is.null(nrow(numEdges))) { 
      exp <- rep(0, length(lists))
    } else {
      exp <- unname(apply(numEdges, 2, sum))
    }
    exp <- c(exp, nrow(edges))
    res[i,] <- exp
  }
  return(res)
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


duplicateEdges <- function(nodes2dup, nodes, edges) {
   for (j in 1:length(nodes2dup)) {
      pos1 <- grep(paste0("^", nodes2dup[j], "$"), edges$from)
      pos2 <- grep(paste0("^", nodes2dup[j], "$"), edges$to)

      len <- length(grep(paste0("^", nodes2dup[j], "$"), nodes$label))
      suf <- seq(1, len, 1)
      
      dup1 <- edges[rep(pos1,len), ]
      dup2 <- edges[rep(pos2,len), ]
      
      if (nrow(dup1)>0) {
        dup1 <- dup1[order(dup1$from, dup1$to), ]
        dup1$from <- paste0(dup1$from, "_", suf)
        dup1$interaction_typeA <- unlist(strsplit(edges$interaction_typeA[pos1], split = "/"))
      }
      if (nrow(dup2)>0) {
        dup2 <- dup2[order(dup2$from, dup2$to), ]
        dup2$to <- paste0(dup2$to, "_", suf)
        dup2$interaction_typeB <- unlist(strsplit(edges$interaction_typeB[pos2], split = "/"))
      }
      edges <- edges[-c(pos1, pos2), ]
      edges <- rbind.data.frame(edges, dup1, dup2)
    }
  return(edges)
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

  
