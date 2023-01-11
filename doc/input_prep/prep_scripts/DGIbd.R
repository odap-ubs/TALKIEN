require(dplyr)
require(tidyr)

load("app/data/geneAnnot.Rdata")

destFile <- "prep_data/dgbiInteractions.tsv"
download.file("https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv", destFile)
dgbi <- read.delim("prep_data/dgbiInteractions.tsv", header = T, stringsAsFactors = F)

idResource$gene.symbol <- NULL
dgbi$entrez_id <- as.character(dgbi$entrez_id)
dgbi$gene_claim_name <- NULL

dgbi <- left_join(dgbi, idResource, by = c("entrez_id" = "entrez.id"))
names(dgbi) <- gsub("\\.", "_", names(dgbi))

dgbi <- dgbi[, c("gene_name", "entrez_id", "uniprot_id", "protein_id", "string_id", "interaction_claim_source", "interaction_types", "drug_claim_name", "drug_claim_primary_name", "drug_name", "drug_concept_id", "interaction_group_score", "PMIDs")]

names(dgbi)[names(dgbi)%in%"uniprot_id"] <- "UniProt_id"

dgbi$drug_claim_primary_name[dgbi$drug_claim_primary_name==""] <- dgbi$drug_claim_name[dgbi$drug_claim_primary_name==""]
dgbi$drug_claim_primary_name <- tolower(dgbi$drug_claim_primary_name)
dgbi$drug_name[dgbi$drug_name==""] <- dgbi$drug_claim_primary_name[dgbi$drug_name==""]
dgbi <- dplyr::select(dgbi, -drug_claim_primary_name, -drug_claim_name) %>% filter(gene_name != "") %>% filter(entrez_id != "") %>% filter(UniProt_id != "") %>% filter(protein_id != "") %>% filter(string_id != "")

dupDrugs <- unique(dgbi$drug_name[duplicated(dgbi$drug_name)])
notdupDrugs <- dgbi$drug_name[!duplicated(dgbi$drug_name)]

dgbi_dupDrugs <- data.frame(matrix(NA, ncol = ncol(dgbi), nrow = nrow(dgbi)))
names(dgbi_dupDrugs) <- names(dgbi)
count <- 1

for (i in 1:length(dupDrugs)) {
  print(dupDrugs[i])
  tmp <- dgbi[dgbi$drug_name%in%dupDrugs[i], ]
  tmp2 <- unique(tmp$gene_name[duplicated(tmp$gene_name, tmp$drug_name)])
  if (length(tmp2)>0) {
    for (j in 1:length(tmp2)) {
      tmp3 <- tmp[grep(tmp2[[j]], tmp$gene_name), ]
      int_type <- paste(unique(tmp3$interaction_types), collapse = ",")
      int_claim <- paste(unique(tmp3$interaction_claim_source), collapse = ",")
      pmids <- paste(unique(unlist(strsplit(tmp3$PMIDs, split = ","))), collapse = ",")
      tmp3 <- tmp3[1, ]
      tmp3$interaction_types <- int_type
      tmp3$interaction_claim_source <- int_claim
      tmp3$PMIDs <- pmids
      dgbi_dupDrugs[count, ] <- tmp3
      rm(tmp3)
      count <- count + 1
    }
  } else {
    dgbi_dupDrugs[count:(count+nrow(tmp)-1), ] <- tmp
    count <- count + nrow(tmp)
  }
  rm(tmp, tmp2)
}
dgbi_dupDrugs <- dgbi_dupDrugs[1:count, ]
dgbi_notdup <- dgbi[dgbi$drug_name%in%notdupDrugs[!notdupDrugs%in%dupDrugs], ]

dgbi_final <- rbind.data.frame(dgbi_dupDrugs, dgbi_notdup)
dgbi_final$interaction_types <- gsub("^,", "", dgbi_final$interaction_types)
dgbi_final$interaction_types <- gsub(",$", "", dgbi_final$interaction_types)
dgbi_final$interaction_claim_source <- gsub("^,", "", dgbi_final$interaction_claim_source)
dgbi_final$interaction_claim_source <- gsub(",$", "", dgbi_final$interaction_claim_source)

dgbi <- dgbi_final %>% arrange(desc(interaction_group_score), drug_name) %>% relocate(drug_name, .after = string_id) %>% relocate(interaction_group_score, .after = drug_name) %>% relocate(interaction_claim_source, .before = PMIDs)

save(dgbi, file = "app/data/DGBidb.Rdata", version = 2)
