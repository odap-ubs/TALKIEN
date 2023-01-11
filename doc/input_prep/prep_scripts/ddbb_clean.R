options(stringsAsFactors = F)
require(igraph)
require(dplyr)

# functions
symbol2string <- function(netlist, resource = idResource) {
  netlist$partner_a <- idResource$string.id[match(netlist$partner_a, idResource$gene.symbol)]
  netlist$partner_b <- idResource$string.id[match(netlist$partner_b, idResource$gene.symbol)]
  return(netlist)
}

igraphGen <- function(net, labs = idResource) {
  graph <- graph_from_data_frame(net) %>% simplify(remove.loops = T)
  V(graph)$labels <- labs$gene.symbol[match(names(V(graph)), labs$string.id)]
  return(graph)
}

# setting dest files
cellChat <- "prep_data/CellChatDB.rda"
cellphone <- "prep_data/cellphoneDB_int.txt"
cellphoneInfo <- "prep_data/cellphoneDB_prots.txt"
icellnet <- "prep_data/icellnetDB.csv"
ramilowski <- "prep_data/ramilowski_int.txt"
ramilowskiLoc <- "prep_data/ramilowski_prots.txt"

download.file("https://github.com/sqjin/CellChat/blob/master/data/CellChatDB.human.rda", cellChat)
download.file("https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/interaction_curated.csv", cellphone)
download.file("https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv", cellphoneInfo)
download.file("https://github.com/soumelis-lab/ICELLNET/blob/master/data/ICELLNETdb.tsv", icellnet)
download.file("https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/PairsLigRec.txt", ramilowski)
download.file("https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/SubcelLoc.Ages.Proteins.txt", ramilowskiLoc)

load("app/data/geneAnnot.Rdata")
load(cellChat)
cellphone <- read.delim(cellphone, header = T, sep = ",")
cellphoneInfo <- read.delim(cellphoneInfo, header = T, sep = ",")
icellnet <- read.delim(icellnet, header = T, sep = ";")
ramilowski <- read.delim(ramilowski, header = T)
ramilowskiLoc <- read.delim(ramilowskiLoc, header = T)

# setting ligand-receptors properly 
# hardcoding HLAA-C to P04439, P01889, P10321
cellphoneInfo <- cellphoneInfo[cellphoneInfo$uniprot%in%idResource$uniprot.id, ]
cellphoneInfo <- cellphoneInfo[!(cellphoneInfo$secreted=="False" & cellphoneInfo$receptor=="False"), ]
cellphoneInfo <- cellphoneInfo[!cellphoneInfo$secreted=="" | !cellphoneInfo$receptor=="", ]
cpnodes <- unique(c(cellphone$partner_a, cellphone$partner_b))
cpnodes <- cpnodes[cpnodes%in%cellphoneInfo$uniprot]

cellphoneInfo <- cellphoneInfo[cellphoneInfo$uniprot%in%cpnodes, ]
cellphoneInfo <- cellphoneInfo[, c("uniprot", "protein_name", "secreted", "receptor")]
cellphoneInfo$uniprot[grep("HLA", cellphoneInfo$uniprot)] <- c("P04439", "P01889", "P10321")

cellphone <- cellphone[cellphone$partner_a%in%cpnodes & cellphone$partner_b%in%cpnodes, ]
cellphone <- cellphone[, c("partner_a", "partner_b")]
cellphone$partner_a[grep("HLAA", cellphone$partner_a)] <- "P04439"
cellphone$partner_b[grep("HLAA", cellphone$partner_a)] <- "P04439"
cellphone$partner_a[grep("HLAB", cellphone$partner_a)] <- "P01889"
cellphone$partner_b[grep("HLAB", cellphone$partner_a)] <- "P01889"
cellphone$partner_a[grep("HLAC", cellphone$partner_a)] <- "P10321"
cellphone$partner_b[grep("HLAC", cellphone$partner_a)] <- "P10321"

cellphone_amb <- cellphoneInfo$uniprot[cellphoneInfo$secreted=="True" & cellphoneInfo$receptor=="True" | cellphoneInfo$receptor==""]
cellphone_clear <- cellphoneInfo[!(cellphoneInfo$uniprot%in%cellphone_amb), ]
cellphone_clear$class <- ifelse(cellphone_clear$secreted=="True", "Ligand", "Receptor")
cellphone_clear$symbol <- idResource$gene.symbol[match(cellphone_clear$uniprot, idResource$uniprot.id)]
cellphone_clear <- cellphone_clear[, c("uniprot", "symbol", "class")]
cellphone_amb <- idResource$gene.symbol[match(cellphone_amb, idResource$uniprot.id)]

cellphone$partner_a <- cellphone_clear$symbol[match(cellphone$partner_a, cellphone_clear$uniprot)]
cellphone$partner_b <- cellphone_clear$symbol[match(cellphone$partner_b, cellphone_clear$uniprot)]
cellphone <- cellphone[!(is.na(cellphone$partner_a) | is.na(cellphone$partner_b)), ]
cellphone$class_a <- cellphone_clear$class[match(cellphone$partner_a, cellphone_clear$symbol)]
cellphone$class_b <- cellphone_clear$class[match(cellphone$partner_b, cellphone_clear$symbol)]
# cellphone cleaned !!!

icellnet1_1 <- icellnet[, c("Ligand.1", "Receptor.1")]
icellnet1_2 <- icellnet[, c("Ligand.1", "Receptor.2")]
icellnet1_3 <- icellnet[, c("Ligand.1", "Receptor.3")]
icellnet2_1 <- icellnet[, c("Ligand.2", "Receptor.1")]
icellnet2_2 <- icellnet[, c("Ligand.2", "Receptor.2")]
icellnet2_3 <- icellnet[, c("Ligand.2", "Receptor.3")]
# achtung, self loops found
names(icellnet1_1) <- names(icellnet1_2) <- names(icellnet1_3) <- names(icellnet2_1) <- names(icellnet2_2) <- names(icellnet2_3) <- c("Ligand", "Receptor")
icellnet <- rbind.data.frame(icellnet1_1, icellnet1_2, icellnet1_3, icellnet2_1, icellnet2_2, icellnet2_3)
icellnet <- icellnet[!(icellnet$Ligand == "" | icellnet$Receptor == ""), ]
icellnet <- icellnet[!(icellnet$Ligand==icellnet$Receptor), ]
icellnet <- icellnet[icellnet$Ligand%in%idResource$gene.symbol & icellnet$Receptor%in%idResource$gene.symbol, ]
rm(icellnet1_1, icellnet1_2, icellnet1_3, icellnet2_1, icellnet2_2, icellnet2_3)

icellnet_amb <- unique(icellnet$Ligand)[unique(icellnet$Ligand)%in%unique(icellnet$Receptor)]
icellnet_clear <- unique(c(icellnet$Ligand, icellnet$Receptor))
icellnet_clear <- data.frame(symbol=icellnet_clear[!icellnet_clear%in%icellnet_amb])
icellnet_clear$class <- ifelse(icellnet_clear$symbol%in%icellnet$Ligand, "Ligand", "Receptor")
# icellnet cleaned !!!

ramilowski_sub <- unique(sort(c(grep("secreted", ramilowskiLoc$Consensus_SL), grep("plasma membrane", ramilowskiLoc$Consensus_SL))))
ramilowskiLoc <- ramilowskiLoc[ramilowski_sub, ]
ramilowskiLoc <- ramilowskiLoc[ramilowskiLoc$ApprovedSymbol%in%idResource$gene.symbol, ]
ramilowski <- ramilowski[ramilowski$Ligand.ApprovedSymbol%in%ramilowskiLoc$ApprovedSymbol | ramilowski$Receptor.ApprovedSymbol %in% ramilowskiLoc$ApprovedSymbol, ]
ramilowski <- ramilowski[!grepl("EXCLUDED", ramilowski$Pair.Evidence), ]
ramilowskiLoc <- ramilowskiLoc[ramilowskiLoc$ApprovedSymbol%in%unique(c(ramilowski$Ligand.ApprovedSymbol, ramilowski$Receptor.ApprovedSymbol)), ]
rm(ramilowski_sub)

ramilowski_amb <- ramilowskiLoc$ApprovedSymbol[grepl("plasma membrane", ramilowskiLoc$Consensus_SL) & grepl("secreted", ramilowskiLoc$Consensus_SL)]
ramilowski_clear <- data.frame(symbol=ramilowskiLoc$ApprovedSymbol[!ramilowskiLoc$ApprovedSymbol%in%ramilowski_amb])
ramilowski_clear$class <- ramilowskiLoc$Consensus_SL[match(ramilowski_clear$symbol, ramilowskiLoc$ApprovedSymbol)]
ramilowski_clear$class <- ifelse(grepl("secreted", ramilowski_clear$class), "Ligand", "Receptor")
# ramilowsky cleaned !!!

# hardcoding some typos detected: activins, inhibins, complexes...
cellchat <- CellChatDB.human$interaction
cellchat_lig <- unique(cellchat$ligand)
cellchat_lig <- gsub("Activin ", "INH", cellchat_lig) # sort of typo; ligand is inhibin not activin
cellchat_lig <- gsub("Inhibin ", "INH", cellchat_lig)
cellchat_lig <- gsub("IL23 complex", "IL23A", cellchat_lig)
cellchat_lig <- gsub("IL27 complex", "IL27A", cellchat_lig)
cellchat_lig[grep("complex", cellchat_lig)] <- gsub(" COMPLEX ", "", toupper(cellchat_lig[grep("complex", cellchat_lig)]))
cellchat_lig[grep("ITGA", cellchat_lig)] <- c("ITGA4", "ITGB1", "ITGA9")
cellchat_lig <- c(cellchat_lig, "ITGB7")
cellchat_lig <- unique(cellchat_lig)
cellchat_lig <- cellchat_lig[cellchat_lig%in%idResource$gene.symbol]

cellchat_rec <- unique(cellchat$receptor)
cellchat_rec <- unlist(strsplit(cellchat_rec, "_"))
cellchat_rec[cellchat_rec=="GP complex"] <- "GP1BA"
cellchat_rec[cellchat_rec=="IL11R complex"] <- "IL11RA"
cellchat_rec <- c(cellchat_rec, "GP1BB", "GP5", "GP9", "IL11ST")
cellchat_rec[grep(":", cellchat_rec)] <- "CD94"
cellchat_rec <- c(cellchat_rec, "NKG2A", "NKG2C", "NKG2E")
cellchat_rec[grep("receptor", cellchat_rec)] <- "CD8A"
cellchat_rec <- c(cellchat_rec, "CD8B")
cellchat_rec[cellchat_rec%in%"R2"] <- "TGFbR2"
cellchat_rec <- toupper(cellchat_rec)
cellchat_rec <- unique(cellchat_rec)
cellchat_rec <- cellchat_rec[cellchat_rec%in%idResource$gene.symbol]

cellchat_amb <- unique(cellchat_lig)[unique(cellchat_lig)%in%unique(cellchat_rec)]
cellchat_clear <- unique(c(cellchat_lig, cellchat_rec))
cellchat_clear <- data.frame(symbol = cellchat_clear[!cellchat_clear%in%cellchat_amb])
cellchat_clear$class <- ifelse(cellchat_clear$symbol%in%cellchat_rec, "Receptor", "Ligand")
# cellchat cleaned!!!

# now time to clean everything without complete annotation and finally recall ambiguous proteins based on comparisons.
unique_ids <- sort(unique(c(cellchat_clear$symbol, cellphone_clear$symbol, icellnet_clear$symbol, ramilowski_clear$symbol)))
ddbbs_annot <- as.data.frame(matrix(nrow=length(unique_ids), ncol = 4))
row.names(ddbbs_annot) <- unique_ids
names(ddbbs_annot) <- c("cellchat", "cellphone", "icellnet", "ramilowsky")
ddbbs_annot$cellchat <- cellchat_clear$class[match(unique_ids, cellchat_clear$symbol)]
ddbbs_annot$cellphone <- cellphone_clear$class[match(unique_ids, cellphone_clear$symbol)]
ddbbs_annot$icellnet <- icellnet_clear$class[match(unique_ids, icellnet_clear$symbol)]
ddbbs_annot$ramilowski <- ramilowski_clear$class[match(unique_ids, ramilowski_clear$symbol)]

# let's look at the discordancies between ddbbs
cons <- unname(unlist(apply(ddbbs_annot, 1, function(class) length(grep("Receptor", class))+length(grep("Ligand", class))*-1)))
ddbbs_annot$consensus <- NA
ddbbs_annot$consensus[cons>0] <- "Receptor"
ddbbs_annot$consensus[cons<0] <- "Ligand"

# some of the ambiguous proteins are not lost because they're found in other ddbbs
table(cellchat_amb%in%row.names(ddbbs_annot)) # 33 recalled
table(cellphone_amb%in%row.names(ddbbs_annot)) # 37 recalled
table(icellnet_amb%in%row.names(ddbbs_annot)) # 3 recalled
table(ramilowski_amb%in%row.names(ddbbs_annot)) # 94 recalled

# writing table
write.table(ddbbs_annot, "suppl/ddbbs_annot.csv", row.names = T, col.names = T, sep = ",", quote = F)


# now time to clean interactions: within ddbb already cleaned, between ddbb need to filter discordant annotations
icellnet <- icellnet[!(icellnet$Ligand%in%icellnet_amb | icellnet$Receptor%in%icellnet_amb), ]
cellphone <- cellphone[!(cellphone$partner_a%in%cellphone_amb | cellphone$partner_b%in%cellphone_amb), ]
ramilowski <- ramilowski[!(ramilowski$Ligand.ApprovedSymbol%in%ramilowski_amb | ramilowski$Receptor.ApprovedSymbol%in%ramilowski_amb), ]
ramilowski <- ramilowski[, c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

cellchat$ligand[grep("Activin", cellchat$ligand)] <- gsub("Activin ", "INH", cellchat$ligand[grep("Activin", cellchat$ligand)])
cellchat$ligand[grep("Inhibin", cellchat$ligand)] <- gsub("Inhibin ", "INH", cellchat$ligand[grep("Inhibin", cellchat$ligand)])
cellchat$ligand[grep("complex", cellchat$ligand)] <- gsub(" complex", "", cellchat$ligand[grep("complex", cellchat$ligand)])
cellchat$ligand <- toupper(cellchat$ligand)
cellchat$ligand <- gsub(" ", "", cellchat$ligand)

cellchat$receptor[grep("_R2$", cellchat$receptor)] <- gsub("R2", "TGFBR2", cellchat$receptor[grep("_R2$", cellchat$receptor)])
cellchat$receptor[grep("GP complex", cellchat$receptor)] <- gsub("GP complex", "GP1BA_GP1BB_GP5_GP9", cellchat$receptor[grep("GP complex", cellchat$receptor)])
cellchat$receptor[grep("IL11R complex", cellchat$receptor)] <- gsub("IL11R complex", "IL11RA_IL6ST", cellchat$receptor[grep("IL11R complex", cellchat$receptor)])
cellchat$receptor <- gsub(":", "_", cellchat$receptor)
cellchat$receptor <- toupper(cellchat$receptor)

ligs <- strsplit(cellchat$ligand, "_")
recs <- strsplit(cellchat$receptor, "_")
ligs <- do.call(rbind.data.frame, ligs)
recs <- do.call(rbind.data.frame, recs)
names(ligs) <- c("ligand1", "ligand2")
names(recs) <- c("receptor1", "receptor2", "receptor3", "receptor4")

cc_int11 <- cbind(ligs$ligand1, recs$receptor1)
cc_int12 <- cbind(ligs$ligand1, recs$receptor2)
cc_int13 <- cbind(ligs$ligand1, recs$receptor3)
cc_int14 <- cbind(ligs$ligand1, recs$receptor4)

cc_int21 <- cbind(ligs$ligand2, recs$receptor1)
cc_int22 <- cbind(ligs$ligand2, recs$receptor2)
cc_int23 <- cbind(ligs$ligand2, recs$receptor3)
cc_int24 <- cbind(ligs$ligand2, recs$receptor4)

cellchat <- rbind(cc_int11, cc_int12, cc_int13, cc_int14, cc_int21, cc_int22, cc_int23, cc_int24)
cellchat <- as.data.frame(cellchat)
names(cellchat) <- c("Ligand", "Receptor")
cellchat <- cellchat[!duplicated(cellchat), ]

names(ramilowski) <- names(icellnet) <- names(cellchat) <- c("partner_a", "partner_b")
ramilowski$class_a <- "Ligand"
ramilowski$class_b <- "Receptor"

icellnet$class_a <- "Ligand"
icellnet$class_b <- "Receptor"

cellchat$class_a <- "Ligand"
cellchat$class_b <- "Receptor"
### prepared individual ligand-receptor networks

ddbbs_merged <- ddbbs_annot[!is.na(ddbbs_annot$consensus), ]
ddbbs_merged <- ddbbs_merged[unname(unlist(apply(ddbbs_merged, 1, function(nodes) length(table(nodes))==1))), ]

icellnet2merge <- icellnet[icellnet$partner_a%in%row.names(ddbbs_merged) & icellnet$partner_b%in%row.names(ddbbs_merged), ]
icellnet2merge <- icellnet2merge[icellnet2merge$partner_a%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Ligand"] & icellnet2merge$partner_b%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Receptor"], ]
icellnet2merge <- symbol2string(icellnet2merge)

ramilowski2merge <- ramilowski[ramilowski$partner_a%in%row.names(ddbbs_merged) & ramilowski$partner_b%in%row.names(ddbbs_merged), ]
ramilowski2merge <- ramilowski2merge[ramilowski2merge$partner_a%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Ligand"] & ramilowski2merge$partner_b%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Receptor"], ]
ramilowski2merge <- symbol2string(ramilowski2merge)

cellchat2merge <- cellchat[cellchat$partner_a%in%row.names(ddbbs_merged) & cellchat$partner_b%in%row.names(ddbbs_merged), ]
cellchat2merge <- cellchat2merge[cellchat2merge$partner_a%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Ligand"] & cellchat2merge$partner_b%in%row.names(ddbbs_merged)[ddbbs_merged$consensus=="Receptor"], ]
cellchat2merge <- symbol2string(cellchat2merge)

cellphone2merge <- cellphone[cellphone$partner_a%in%row.names(ddbbs_merged) & cellphone$partner_b%in%row.names(ddbbs_merged), ]
cpRec <- unique(c(cellphone2merge$partner_a[cellphone2merge$class_a=="Receptor"], cellphone2merge$partner_b[cellphone2merge$class_b=="Receptor"]))
cpLig <- unique(c(cellphone2merge$partner_a[cellphone2merge$class_a=="Ligand"], cellphone2merge$partner_b[cellphone2merge$class_b=="Ligand"]))
cpTot <- c(cpRec, cpLig)
cellphone2merge <- cellphone2merge[cellphone2merge$partner_a%in%cpTot & cellphone2merge$partner_b%in%cpTot, ]
cellphone2merge <- symbol2string(cellphone2merge)

ddbbs_merged$vertices <- idResource$string.id[match(row.names(ddbbs_merged), idResource$gene.symbol)]
ddbbs_merged <- ddbbs_merged[, c("vertices", "consensus")]
names(ddbbs_merged)[2] <- "class"

crossTalkNet <- rbind.data.frame(cellchat2merge, cellphone2merge, icellnet2merge, ramilowski2merge)
crossTalkNet <- crossTalkNet[!duplicated(crossTalkNet), ]
crossTalknet <- symbol2string(crossTalkNet)

crossTalkNetIgraph <- igraphGen(crossTalkNet)
crossTalkNetCellChat <- igraphGen(cellchat2merge)
crossTalkNetCellPhone <- igraphGen(cellphone2merge)
crossTalkNetiCellNet <- igraphGen(icellnet2merge)
crossTalkNetRam <- igraphGen(ramilowski2merge)


save(cellchat, cellphone, icellnet, ramilowski, ddbbs_annot, ddbbs_merged, file = "app/data/LigRecInts.Rdata", version = 2)
save(crossTalkNetIgraph, crossTalkNetCellChat, crossTalkNetCellPhone, crossTalkNetiCellNet, crossTalkNetRam, file = "app/data/ctnets.Rdata", version = 2)
