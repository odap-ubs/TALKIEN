expInts <- function(dat, graph, lists) {
  
  start <- Sys.time()
  res <- data.frame(matrix(NA, nrow = 1000, ncol = length(lists)+1))
  names(res) <- c(lists, "total")
  nCond <- unique(dat[,2])
  for (i in 1:1000) {
    nodeList <- list()
    for (j in 1:length(nCond)) {
      nodeList[[j]] <- sample(idResource$gene.symbol, replace = F, size = table(dat[,2])[[j]])
    }
    names(nodeList) <- nCond
    nodeList <- lapply(nodeList, function(cond) { 
      ind <- cond%in%row.names(ddbbs_merged)
      id <- cond[ind]
      class <- ddbbs_merged$class[match(id, row.names(ddbbs_merged))]
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
    nodes$stringId <- idResource$string.id[match(nodes$id, idResource$gene.symbol)]
    nodeIds <- names(V(graph))[names(V(graph))%in%nodes$stringId]
    nodes <- nodes[nodes$stringId%in%nodeIds, ]
    net <- induced_subgraph(as.undirected(graph), nodeIds) %>%
    net <- set.vertex.attribute("class", value = nodes$class[match(names(V(net)), nodes$stringId)]) %>%
    net <- set.vertex.attribute("list", value = nodes$list[match(names(V(net)), nodes$stringId)])
    edges <- as.data.frame(as_edgelist(net))
    from <- nodes[match(edges[,1], nodes$stringId), c("list", "class")]
    to <- nodes[match(edges[,2], nodes$stringId), c("list", "class")]
    edges <- cbind.data.frame(edges, from, to)
    names(edges) <- c("from", "to", "list_from", "class_from", "list_to", "class_to")
    #edges_crt <- filter(edges, list_from !=list_to & class_from != class_to)
    
    numEdges <- sapply(lists, function(x) paste(edges$list_from, edges$list_to, sep = " X ") %in% x | paste(edges$list_to, edges$list_from, sep = " X ") %in% x)
    res[i,1:(ncol(res)-1)] <- unname(apply(numEdges, 2, sum))
    res$total[i] <- nrow(numEdges)
  }
  
  #pval1=as.character(sum(abs(res1) >= iEdges)/(1000))
  #pval2=as.character(sum(abs(res2) >= iEdges_crt)/(1000))
  #if (pval1 == 0) pval1 <- "< 0.001"
  #if (pval2 == 0) pval2 <- "< 0.001"
  print(Sys.time()-start)
  return(res)
}
