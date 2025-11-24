# Function: Differentially expressed circRNAs and pathway enrichment.
# Programmer: Zhangli Yuan, Yi Shi, 20251124.

##### DESeq2 ######
library(DESeq2)
library(tidyverse)
library(ggplot2)

counts <- read.csv("GSE158596_counts.csv")
circRNAID <- counts$circRNA_ID
counts <- counts[, -1]
rownames(counts) <- circRNAID

condition <- factor(c(rep('notcancer', 16), rep('cancer', 78)),
                    levels = c('notcancer', 'cancer'))
colData <- data.frame(row.names = colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "notcancer")
dds_2 <- DESeq(dds)
res <- results(dds_2)
save(res,file = "DEG.Rda")

DEG <- as.data.frame(res)    
exp <- counts

logFC_cutoff <- 0.15
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)


library(pheatmap)
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
write.csv(c, "Diff_circRNA_id.csv", row.names = TRUE)
exp_diff0 <- cbind(circRNA_ID = rownames(exp_diff), exp_diff)
write.csv(exp_diff0, "Diff_counts.csv", row.names = FALSE)
write.csv(a,file = "cancer_up.csv")
write.csv(b,file = "cancer_down.csv")

#------------
annotation_col <- colData
a <- filter(annotation_col,condition == 'cancer')
b <- filter(annotation_col,condition == 'notcancer')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)

pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,              
         show_colnames =F,                
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,              
         cluster_rows = T,             
         fontsize = 10,           
         fontsize_row=3,
         fontsize_col=3)
dev.off()

#-----------
summary(res)
res <- res[order(res$padj), ]
circRNA <- rownames(res)

resdata <- as_tibble(res) %>%
  mutate(sig = case_when(
    log2FoldChange >= 0.15 & padj < 0.05 ~ 'up(p.adj < 0.05, log2FC >= 0.15)',
    log2FoldChange <= -0.15 & padj < 0.05 ~ 'down(p.adj < 0.05, log2FC <= -0.15)',
    TRUE ~ 'no diff'
  ))
resdata <- as_tibble(res) %>%
  mutate(sig = case_when(
    log2FoldChange >= 0.15 & padj < 0.05 ~ 'up',
    log2FoldChange <= -0.15 & padj < 0.05 ~ 'down',
    TRUE ~ 'normal'
  ))

resdata <- cbind(circRNA, resdata)
write.csv(resdata, 'resdata.csv', row.names = FALSE)

filtered_data <- subset(resdata, abs(-log(padj, 10)) < 25)

ggplot(data = filtered_data, mapping = aes(x = log2FoldChange, y = -log(padj, 10))) +
  geom_point(mapping = aes(color = sig), alpha = 0.6, size = 1) +
  scale_color_manual(values = c('blue2', 'gray30', 'red2')) +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = c(0.85, 0.8),
        legend.title = element_blank(),
        legend.key = element_rect(fill = 'transparent'),
        legend.box.background = element_rect(color = 'gray', size = 0.5),
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.15, 0.15), color = 'gray', size = 01, lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 01, lwd = 0.5, lty = 2) +
  labs(x = 'log2FoldChange', y = '-log10(p.value)', color = NA) +
  ggtitle('GSE158596')

ggsave('GSE158596_volcano.png')
dev.off()



##### Pathway Enrichment ######
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

DEG <- as.data.frame(res)

DEG <- DEG[circID$circRNA_ID,]
DEG <- rownames_to_column(DEG,"circRNA")

geneid <- read.csv("gene_id_trans.csv")

duplicated_rows <- duplicated(geneid$circRNA_ID)
gene_id_uni <- geneid[!duplicated_rows, ]

row.names(gene_id_uni) <- gene_id_uni$circRNA_ID
gene_deg <- gene_id_uni[DEG$circRNA,]
b <- is.na(gene_deg$circRNA_ID)
DEG <- DEG[!b,]
gene_deg <- gene_deg[!b,]
write.csv(gene_deg,"Diff_gene_id.csv", row.names = FALSE)

row.names(gene_deg) <- NULL
DEG <- cbind(DEG, gene_deg)
dupli_rows <- duplicated(DEG$SYMBOL)
DEG <- DEG[!dupli_rows, ]

genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL") 

# GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_DEG_final.Rda")
write.csv(ego_res,"GO_DEG_result.csv")

#------------
go_gene <- ego_res$geneID
go_gene_list <- NULL
for (i in go_gene){
  split <- strsplit(i, "/")
  for (j in split){
    go_gene_list <- c(go_gene_list, j)
  }
  
}
go_genes <- data.frame(go_gene_list)
dupli_rows <- duplicated(go_genes$go_gene_list)
go_genes <- go_genes[!dupli_rows, ]
go_genes <- data.frame(go_genes)
write.csv(go_genes, "go_genes_symbol.csv",row.names = FALSE)

library(tidyverse)
library(conflicted)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

geneid <- read.csv("Diff_gene_id.csv")
go_gene <- go_genes
go_gene <- read.csv("go_genes_symbol.csv")
gene_id_go <- geneid %>% filter(SYMBOL %in% go_gene$SYMBOL)
dupli_rows <- duplicated(gene_id_go$circRNA_ID)
gene_id_go <- gene_id_go[!dupli_rows, ]
write.csv(gene_id_go,"go_gene_id.csv",row.names = FALSE) 

counts <- read.csv("GSE158596_counts.csv")
go_counts <- counts %>% filter(circRNA_ID %in% gene_id_go$circRNA_ID)
write.csv(go_counts,"go_counts.csv",row.names = FALSE)

#----------
library(tidyverse)
library(dplyr)
count <- read.csv("GSE158596_counts.csv", row.names= 1)

cancer_State <- NULL
for (i in 1:16){
  cancer_State <- c(cancer_State,0)
}
for (i in 1:78){
  cancer_State <- c(cancer_State,1)
}

library(DescTools)
gene_correlations <- apply(count, 1, function(gene_expr) {
  cor(gene_expr, cancer_State, method = "spearman")
})
correlation <- data.frame(gene_correlations)
write.csv(correlation, "all_circRNA_aml_correlation.csv",row.names = TRUE)


go_count <- read.csv("go_counts.csv", row.names= 1)
diff_count <- read.csv("Diff_counts.csv", row.names= 1)

gene_correlations <- apply(diff_count, 1, function(gene_expr) {
  cor(gene_expr, cancer_State, method = "spearman")
})
correlation <- data.frame(gene_correlations)
write.csv(correlation, "circRNA_aml_correlation.csv",row.names = TRUE)
correlation_abs <- correlation %>% mutate_all(abs)
write.csv(correlation_abs, "circRNA_aml_correlation_abs.csv",row.names = TRUE)

correlation_abs <- rownames_to_column(correlation_abs)
correlation_abs <- correlation_abs[order(correlation_abs$gene_correlations, decreasing = TRUE), ]
row.names(correlation_abs) <- NULL
write.csv(correlation_abs, "diff_circRNA_aml_correlation_abs.csv",row.names = FALSE)

#---------
library(readr)
go_gene <- ego_res$geneID
max_comma_count <- max(str_count(ego_res$geneID, "/"))
cat("The maximum comma count:", max_comma_count, "\n")
library(tidyr)
geneID <- separate(ego_res, geneID, into = paste0("gene_", 1:33), sep = "/")
write.csv(geneID,"go_p_symbol.csv",row.names = TRUE)

gene_id_trans <- read.csv("go_gene_id.csv")
library(readxl)
library(tidyverse)
library(dplyr)
correlation <- read.csv("circRNA_aml_correlation_abs.csv")
correlation <- correlation[order(correlation$circRNA_ID), ]
gene_id_trans <- gene_id_trans[order(gene_id_trans$circRNA_ID), ]
row.names(correlation) <- NULL
row.names(gene_id_trans) <- NULL
gene_correlation <- merge(correlation,gene_id_trans,by = "circRNA_ID")
write.csv(gene_correlation, "go_gene_correlation.csv", row.names = FALSE)

go_symbol <- read_excel("go_symbol.xlsx")
go_1 <- t(go_symbol[1,])
go_1_id <- gene_correlation %>% filter(SYMBOL %in% go_1)
go_1_id <- go_1_id[order(go_1_id$gene_correlations, decreasing = TRUE), ]
row.names(go_1_id) <- NULL

go_id <- go_1_id
group_col <- rep(1, nrow(go_id))
go_id$go_path <- group_col

go_circRNA <- go_1_id$circRNA_ID

for (i in 2:70){
  go_1 <- t(go_symbol[i,])
  go_1_id <- gene_correlation %>% filter(SYMBOL %in% go_1) %>% filter(!circRNA_ID %in% go_circRNA)
  go_1_id <- go_1_id[order(go_1_id$gene_correlations, decreasing = TRUE), ]
  row.names(go_1_id) <- NULL
  group_col <- rep(i, nrow(go_1_id))
  go_1_id$go_path <- group_col
  go_id <- rbind(go_id,go_1_id)
  go_circRNA <- c(go_circRNA,go_1_id$circRNA_ID)
}
write.csv(go_id,"go_path_id.csv",row.names = FALSE)


#-------
library(readxl)
library(tidyverse)
library(dplyr)

go_path_id <- read.csv("go_path_id.csv")
path <- unique(go_path_id$go_path)
go_path_10 <- data.frame()

for (i in path){
  go_path_i <- go_path_id[go_path_id$go_path == i,]
  if (nrow(go_path_i) > 5){
    go_path_i <- go_path_i[1:5,]
  }
  go_path_10 <- rbind(go_path_10,go_path_i)
}
go_id <- go_path_10

counts <- read.csv("GSE158596_counts.csv")
go_path_counts <- counts %>% filter(circRNA_ID %in% go_id$circRNA_ID)
write.csv(go_path_counts,"go_path_counts.csv",row.names = FALSE)

#------
library(readxl)
library(tidyverse)
library(dplyr)

go <- read.csv("GO_DEG_result.csv",row.names = 1)
table(go[,1])

paths <- c(1,2,3,4,5,6,8,9,10,12,16,20,21,22,23,24,26,30,33,39,41,45,46,47,49,50,53,54,57,62,63,64,65,68)
go <- go[paths,]

table(go[,1])
go_MF<-go[go$ONTOLOGY=="MF",]
go_CC<-go[go$ONTOLOGY=="CC",]
go_BP<-go[go$ONTOLOGY=="BP",]
go_enrich_df<-data.frame(ID=c(go_BP$ID, go_CC$ID, go_MF$ID),
                         Description=c(go_BP$Description, go_CC$Description, go_MF$Description),
                         GeneNumber=c(go_BP$Count, go_CC$Count, go_MF$Count),
                         type=factor(c(rep("biological process", 21), rep("cellular component", 3),rep("molecular function",10)),levels=c("molecular function", "cellular component", "biological process")))

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels <- go_enrich_df$Description
names(labels) = rev(1:nrow(go_enrich_df))

CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
p
dev.off()
