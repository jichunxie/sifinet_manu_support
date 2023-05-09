#setwd("~/Desktop/SiFINeT/Result/Experimental/Monoclonal/cellcycle")

annotation <- read.delim("out.txt", header = F, sep = ";")
annotation <- annotation[,c(1,4,5)]
annotation[,2] <- sub('...', '', annotation[,2])
colnames(annotation) <- c("uniprot_id", "annotation", "source")


id <- read.table("uniid-symbol.txt", header = T, sep = "\t")
colnames(id) <- c("uniprot_id", "name")

data <- merge(annotation, id, by = "uniprot_id", all.x = T)

#write.table(data, "annotated_data.txt", quote = F, row.names = F, col.names = T, sep = ";")

rm(annotation, id)

key_annotation <- c("mitotic",
                    #"cell cycle DNA replication",
                    "centriole",
                    "centrosome",
                    "centromere",
                    # "chromosome attachment to the nuclear envelope",
                    "spindle",
                    "chromosome separation",
                    #"cytokine",
                    "cohesin",
                    "cell cycle",
                    "DNA replication",
                    "DNA repair",
                    "microtubule",
                    "nuclear chromosome segregation",
                    #"regulation of cell cycle",
                    "cell division",
                    "cell growth")

out <- c()

for (i in 1:nrow(data)){
  out[i] <- sum(sapply(key_annotation, function(x){grepl(x, data$annotation[i], fixed=TRUE)}))
}
feature_genes <- unique(data$name[out>0])

gene_db <- readRDS("gene_db.rds")
gene_st <- readRDS("gene_st.rds")
all_cc_gene <- union(gene_db, gene_st)

# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  return(humanx)
}

mouse.cc.genes <- convertHumanGeneList(all_cc_gene)

joined_cc_gene <- intersect(tolower(feature_genes), tolower(mouse.cc.genes))

saveRDS(joined_cc_gene, "ccgene.rds")
