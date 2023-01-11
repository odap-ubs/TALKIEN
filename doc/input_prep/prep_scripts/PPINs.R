require(dplyr)
require(tidyr)

load("app/data/geneAnnot.Rdata")

destFile <- "prep_data/9606.protein.links.detailed.v11.5.txt.gz"

download.file("https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz", destFile)

string <- data.table::fread("prep_data/9606.protein.links.detailed.v11.5.txt.gz", header = T, stringsAsFactors = F)

string_hiscore <- dplyr::filter(string, combined_score>700)
string_exp <- dplyr::filter(string, experimental>0)
string_ddbb <- dplyr::filter(string, database>0)

png("suppl/histCurated.png", height = 700, width = 900)
hist(string_exp$experimental, 
     breaks = 20, 
     col = rgb(1,0,0,0.5), 
     main = "Histogram STRING scores", 
     xlab = "experimental score / database score",
     ylab = "freq")

hist(string_ddbb$database, 
     breaks = 20, 
     col = rgb(0,0,1,0.5),
     add = T)

abline(v=300)
legend("topright", 
       legend=c("Experimental score","Database score"), 
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), 
       pt.cex=2, 
       pch=15)
dev.off()

string_cur <- dplyr::filter(string, experimental>300 | database>300)

png("suppl/histCurated_combinedScore.png", height = 700, width = 900)
hist(string_cur$combined_score, 
     breaks = 20, 
     col = "lightgray", 
     main = "Histogram STRING curated", 
     xlab = "combined score",
     ylab = "freq")
dev.off()

string_cur <- dplyr::filter(string_cur, protein1 %in% idResource$string.id & protein2 %in% idResource$string.id)
string_hiscore <- dplyr::filter(string_hiscore, protein1 %in% idResource$string.id & protein2 %in% idResource$string.id)
names(string_cur)[10] <- names(string_hiscore)[10] <- "weight"

ppi_exp <- igraph::graph_from_data_frame(string_cur[,c(1:2,10)])
ppi_hiscore <- igraph::graph_from_data_frame(string_hiscore[,c(1:2,10)])

save(ppi_exp, ppi_hiscore, file = "app/data/PPIs.Rdata", version = 2)
