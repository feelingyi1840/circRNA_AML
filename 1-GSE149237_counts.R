# Function: Generate Counts Matrix Based on CIRI2 Results for GSE149237.
# Programmer: Zhangli Yuan, Yi Shi, 20251124.

library(tidyverse)
library(ggplot2)
library(conflicted)
library(readxl)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

dir_path <- "GSE149237"
out_files <- list.files(dir_path, pattern = ".xlsx$")

file_contents <- list()
for (i in 1:length(out_files)) {
  file_contents[[i]] <- read_excel(out_files[i])
}

circRNA_ID <- list()
for (i in 1:length(out_files)) {
  circRNA_ID[[i]] <- file_contents[[i]]$'circRNA_ID'
}

tmp <- NA
for (i in 1:length(out_files)) {
  tmp <- c(tmp, circRNA_ID[[i]])
}
all_ID <- sort(unique(tmp)[-1])

results <- matrix(0, nrow = length(all_ID), ncol = 26)
for (i in 1:26) {
  counts <- sapply(all_ID, function(x){sum(circRNA_ID[[i]] == x)})
  results[, i] <- counts
}

results_2 <- as_tibble(cbind(results, rowSums(results), index = 1:length(all_ID))) %>%
  filter(V27 >= 1) %>%
  dplyr::select(index) 

all_ID_index <- as_tibble(cbind(all_ID, index = 1:length(all_ID)))
all_ID_index$index <- as.double(all_ID_index$index)
selected_ID <- inner_join(all_ID_index, results_2, by = 'index')
save(selected_ID, file = 'Selected_ID.RData')
write.csv(selected_ID, out_file, row.names = FALSE)

#----
library(tidyverse)
library(ggplot2)
library(conflicted)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')
load(file = 'Selected_ID.RData')

##### gene_id
cancer_1 <- read_excel('cancer_4.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_2 <- read_excel('cancer_4.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_3 <- read_excel('cancer_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_4 <- read_excel('cancer_5.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_5 <- read_excel('cancer_6.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_6 <- read_excel('cancer_6.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_7 <- read_excel('cancer_7.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_8 <- read_excel('cancer_7.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_9 <- read_excel('cancer_8.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_10 <- read_excel('cancer_8.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_11 <- read_excel('cancer_9.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_12 <- read_excel('cancer_9.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_13 <- read_excel('cancer_10.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_14 <- read_excel('cancer_10.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_15 <- read_excel('cancer_11.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_16 <- read_excel('cancer_11.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

notcancer_1 <- read_excel('healthy_4.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_2 <- read_excel('healthy_4.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_3 <- read_excel('healthy_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_4 <- read_excel('healthy_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_5 <- read_excel('healthy_6.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_6 <- read_excel('healthy_6.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_7 <- read_excel('healthy_7.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_8 <- read_excel('healthy_7.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_9 <- read_excel('healthy_8.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_10 <- read_excel('healthy_8.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

notcancer <- bind_rows(notcancer_1, notcancer_2, notcancer_3, notcancer_4, notcancer_5, notcancer_6, notcancer_7, notcancer_8, 
                       notcancer_9, notcancer_10,
                       .id = 'notcancer')
cancer <- bind_rows(cancer_1, cancer_2, cancer_3, cancer_4, cancer_5, cancer_6, 
                    cancer_7, cancer_8, cancer_9, cancer_10,cancer_11, cancer_12, cancer_13, cancer_14, cancer_15, cancer_16,
                    .id = 'cancer')
data <- bind_rows(notcancer, cancer)
write.csv(data, 'GSE149237_data.csv', row.names = FALSE)

#----
gene_id <- bind_rows(notcancer, cancer) %>%
  dplyr::select(circRNA_ID, gene_id) %>%
  dplyr::filter(gene_id != 'n/a')

max_comma_count <- max(str_count(gene_id$gene_id, ","))
cat("The maximum comma count:", max_comma_count, "\n")
library(tidyr)
gene_id <- separate(gene_id, gene_id, into = paste0("gene_", 22), sep = ",")
gene_id_long <- pivot_longer(gene_id, cols = starts_with("gene_"), names_to = "gene_column", values_to = "gene_id") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  dplyr::filter(gene_id != '') %>%
  select(circRNA_ID, gene_id) %>%
  distinct()

write.csv(gene_id_long, 'GSE149237_gene_id.csv', row.names = FALSE)
#----

# merge circRNA data
setwd('D:/prp/test2')
library(tidyverse)
library(ggplot2)
library(conflicted)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')
load(file = 'Selected_ID.RData')

cancer_1 <- read_excel('cancer_4.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_2 <- read_excel('cancer_4.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_3 <- read_excel('cancer_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_4 <- read_excel('cancer_5.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_5 <- read_excel('cancer_6.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_6 <- read_excel('cancer_6.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_7 <- read_excel('cancer_7.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_8 <- read_excel('cancer_7.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_9 <- read_excel('cancer_8.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_10 <- read_excel('cancer_8.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_11 <- read_excel('cancer_9.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_12 <- read_excel('cancer_9.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_13 <- read_excel('cancer_10.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_14 <- read_excel('cancer_10.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_15 <- read_excel('cancer_11.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_16 <- read_excel('cancer_11.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

notcancer_1 <- read_excel('healthy_4.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_2 <- read_excel('healthy_4.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_3 <- read_excel('healthy_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_4 <- read_excel('healthy_5.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_5 <- read_excel('healthy_6.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_6 <- read_excel('healthy_6.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_7 <- read_excel('healthy_7.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_8 <- read_excel('healthy_7.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_9 <- read_excel('healthy_8.1.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_10 <- read_excel('healthy_8.2.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

notcancer <- bind_rows(notcancer_1, notcancer_2, notcancer_3, notcancer_4, notcancer_5, notcancer_6, notcancer_7, notcancer_8, 
                       notcancer_9, notcancer_10,
                       .id = 'notcancer') %>%
  pivot_wider(names_from = notcancer, values_from = `#junction_reads`)
notcancer[is.na(notcancer)] <- 0
colnames(notcancer) <- c('circRNA_ID', paste0('notcancer', '_', 17:26))

cancer <- bind_rows(cancer_1, cancer_2, cancer_3, cancer_4, cancer_5, cancer_6, 
                    cancer_7, cancer_8, cancer_9, cancer_10,cancer_11, cancer_12, cancer_13, cancer_14, cancer_15, cancer_16,
                    .id = 'cancer') %>%
  pivot_wider(names_from = cancer, values_from = `#junction_reads`)
cancer[is.na(cancer)] <- 0
colnames(cancer) <- c('circRNA_ID', paste0('cancer', '_', 79:94))

counts <- full_join(notcancer, cancer, by = 'circRNA_ID')
counts[is.na(counts)] <- 0
counts <- as.data.frame(counts)
write.csv(counts, "GSE149237_counts.csv", row.names = FALSE)
