library(cogena)
library(stringr)

httr::set_config(httr::config(ssl_verifypeer = FALSE))

set.seed(1)
genesetlist <- gmt2list("cellcycle/c2.all.v2023.1.Hs.symbols.gmt")

# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org", mirror = "uswest")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org", mirror = "uswest")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}


path1 <- str_to_title(unlist(genesetlist["KEGG_OXIDATIVE_PHOSPHORYLATION"]))
path2 <- str_to_title(unlist(genesetlist["KEGG_CITRATE_CYCLE_TCA_CYCLE"]))
path3 <- str_to_title(unlist(genesetlist["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]))
path <- c(path1, path2, path3)
path_conv <- convertHumanGeneList(path)
saveRDS(path_conv, "cellcycle/pathgene.rds")

hub1 <- c("AHSG", "SERPINC1", "FGA", "F2", "CP", 
          "ITIH2", "APOA2", "HPX", "PLG", "HRG",
          "TIMP1", "CXCL1", "COL1A2", "MMP1", "AURKA", 
          "UBE2C", "CXCL12", "TOP2A", "ALDH1A1", "PRKACB")
hub2 <- c("ADGRB3", "CCNF", "CKAP2L", "DIAPH3", "OSBPL3", "RERGL")
hub3 <- c("CDK1", "CCNB1", "CCNA2", "AURKA", "CDC20", 
         "AURKB", "TPX2", "BUB1", "CDC45", "MAD2L1",
         "KIF2C", "NCAPG", "DLGAP5", "FOXM1", "CENPF",
         "CENPE", "BUB1B", "TTK", "ASPM", "KIF20A")
hub <- c(hub1, hub2, hub3)
hub_conv <- convertHumanGeneList(hub)
saveRDS(hub_conv, "cellcycle/hubgene.rds")
