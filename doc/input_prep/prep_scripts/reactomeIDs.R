# getting homo sapiens reactome pathways

ezID2path <- as.list(reactome.db::reactomeEXTID2PATHID)
ezID2path <- lapply(ezID2path, function(ids) ids[grepl("R-HSA", ids)])
ezID2path <- ezID2path[lapply(ezID2path,length)>0]

path2name <- as.list(reactome.db::reactomePATHID2NAME)
path2name <- path2name[names(path2name)%in%unname(unlist(ezID2path))]

save(ezID2path, path2name, file = "app/data/reactomeIDs.Rdata", version = 2)
