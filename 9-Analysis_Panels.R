# Function: Analysis of different circRNA panels -- Part 2
# Programmer: Zhangli Yuan, Yi Shi, 20251124.

#### Pathway Enrichment ####
library(org.Hs.eg.db)
library(clusterProfiler)

bin_sign <- read.csv("bin_sign.csv")
geneid <- read.csv("gene_id_trans.csv")
go_res <- data.frame()

for (cbin in 0:17){
  circRNA_bin <- bin_sign %>% filter(cobin == cbin)
  bin_gene <- geneid %>% filter(circRNA_ID %in% circRNA_bin$circRNA_ID) 
  
  duplicated_rows <- duplicated(bin_gene$SYMBOL)
  bin_gene_uni <- bin_gene[!duplicated_rows, ]
  bin_genelist <- bitr(bin_gene_uni$SYMBOL, fromType="SYMBOL",
                       toType="ENTREZID", OrgDb='org.Hs.eg.db')
  
  ego <- enrichGO(gene = bin_genelist$ENTREZID,
                  OrgDb = org.Hs.eg.db,    
                  ont = "all",         
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff =0.05, 
                  qvalueCutoff =0.05,
                  readable = TRUE)
  ego_res <- ego@result
  
  ego_res$bin <- cbin
  go_res <- rbind(go_res,ego_res)
}

write.csv(go_res,"go_bin_result.csv",row.names = FALSE)


#### Venn diagram ####
library(tidyverse)
library(ggplot2)
library(tidyr)
library(conflicted)
library(readxl)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

geneid <- read.csv("gene_id_trans.csv")

panel_topn <- read.csv("control_counts_train+test.csv")$circRNA_ID
panel_pathway <- read.csv("go_path_counts_train+test.csv")$circRNA_ID
panel_3dcluster <- read.csv("circ_3d_counts_train+test.csv")$circRNA_ID
panel_3dradius5 <- read.csv("circ_radius_5_counts_train+test.csv")$circRNA_ID

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(List1 = panel_topn, List2 = panel_pathway, List3 = panel_3dcluster, List4 = panel_3dradius5),
  category.names = c("Panel-TopN", "Panel-Pathway", "Panel-3DG-Cluster", "Panel-3DG-Radius5"),
  filename = "venn_diagram.png",
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  compression = "lzw"
)