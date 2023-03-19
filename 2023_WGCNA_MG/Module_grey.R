setwd("/home/mguidry/repos/Larval-Oyster-CASE-RNA/2023_WGCNA_MG")

#read in module grey genes
grey <- read.csv("module_grey_genes.csv")
dim(grey)
grey.df <- as.data.frame(grey)
names(grey.df) <- c("X", "geneID")
genes <- grey.df$geneID
head(grey.df)
View(grey.df)

# read in Go term file
GOtable <- read.table("STR.GO.list.unknowns.rmsemicolon")
names(GOtable) <- c("geneID", "GO-Term")


merged_df <- merge(grey.df, GOtable, by = "geneID")
grey_GOterms <- data.frame(geneID=merged_df$geneID, GOterm=merged_df$'GO-Term')
View(grey_GOterms)