options(stringsAsFactors = F)

require(data.table)

destFile <- "prep_data/9606.protein.info.v11.5.txt.gz"

download.file("https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz", destFile)

string <- data.table::fread("prep_data/9606.protein.info.v11.5.txt.gz", header = T, stringsAsFactors = F)
uniprot <- read.delim("prep_data/uniprot-download_true_fields_accession_2Cid_2Cprotein_name_2Cgene_na-2022.06.29-09.58.30.80.tsv", header = T)

string <- string[, 1:2]
names(string) <- c("string.id", "gene.symbol")
names(uniprot) <- c("uniprot.id", "protein.id", "protein.name", "gene.symbol")
uniprot$gene.symbol <- gsub(" .*", "", uniprot$gene.symbol)

require(org.Hs.eg.db)
symbol <- unique(uniprot$gene.symbol)
hs <- org.Hs.eg.db
entrez <- select(hs, 
                 keys = symbol,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL")
names(entrez) <- c("gene.symbol", "entrez.id")

idResource <- dplyr::full_join(uniprot, entrez)
idResource <- dplyr::inner_join(idResource, string)
idResource <- idResource[!duplicated(idResource$uniprot.id) & !duplicated(idResource$protein.id) & !duplicated(idResource$gene.symbol) & !duplicated(idResource$entrez.id), ]
idResource <- idResource[complete.cases(idResource), ]

save(idResource, file = "app/data/geneAnnot.Rdata", version = 2)
