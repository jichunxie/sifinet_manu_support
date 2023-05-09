library(cogena)

cc_database <- c("REACTOME_CELL_CYCLE", 
                 "BIOCARTA_CELLCYCLE_PATHWAY",
                 "KEGG_CELL_CYCLE")


data <- gmt2list("c2.all.v2023.1.Hs.symbols.gmt")
a <- names(data)
b1 <- grepl("CELL_CYCLE", a, fixed = TRUE)
b2 <- grepl("CELLCYCLE", a, fixed = TRUE)
b <- b1 | b2
d <- data[b]
db <- names(d) %in% cc_database

db_genes <- unique(unlist(d[db]))
st_genes <- unlist(d[!db])

e <- table(st_genes)
selected_st_genes <- names(e)[e >= 3]

saveRDS(db_genes, "gene_db.rds")
saveRDS(selected_st_genes, "gene_st.rds")