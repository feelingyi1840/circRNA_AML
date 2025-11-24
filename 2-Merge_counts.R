###### Merge GSE158596 and GSE149237 ######

library(tidyverse)

count1 <- read.csv("GSE158596_counts.csv")
count2 <- read.csv("GSE149237_counts.csv")

common_circRNA <- intersect(count1$circRNA_ID, count2$circRNA_ID)
circ_only_1 <- (count1 %>% filter(!circRNA_ID %in% common_circRNA))$circRNA_ID
circ_only_2 <- (count2 %>% filter(!circRNA_ID %in% common_circRNA))$circRNA_ID

col_names <- names(count2)
C_matrix <- matrix(0, nrow = length(circ_only_1), ncol = ncol(count2), dimnames = list(circ_only_1, col_names))
C2 <- as.data.frame(C_matrix)
C2$circRNA_ID <- row.names(C2)
row.names(C2) <- NULL
count2_new <- rbind(count2,C2)                                   
count2_new <- count2_new %>% filter(!circRNA_ID %in% circ_only_2)

count1 <- count1[order(count1$circRNA_ID), ]
count2_new <- count2_new[order(count2_new$circRNA_ID), ]
row.names(count1) <- count1$circRNA_ID
row.names(count2_new) <- count2_new$circRNA_ID
count2_new <- count2_new[,-1]
identical(row.names(count1),row.names(count2_new))
counts <- cbind(count1,count2_new)
row.names(counts) <- NULL

rownames(counts) = counts[,1]
counts = counts[-1]
counts_z <- apply(counts,2,scale)
z_counts <- as.data.frame(counts_z)

write.csv(z_counts,"counts_train+test.csv",row.names = TRUE) 