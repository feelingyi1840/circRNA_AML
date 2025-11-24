###### circRNA Panel Selection ######

library(tidyverse)
library(ggplot2)
library(tidyr)
library(conflicted)
library(readxl)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

N <- 51 

#### Panel-TopN ####
correlation <- read.csv("diff_circRNA_aml_correlation_abs.csv")
control <- correlation[1:N,]
count <- read.csv("counts_cv.csv")
control_counts <- count %>% filter(circRNA_ID %in% control$circRNA_ID)
write.csv(control_counts,"control_3d_counts.csv",row.names = FALSE)


#### Panel-Pathway ####
go_path_id <- read.csv("go_path_id.csv")
path <- unique(go_path_id$go_path)
go_path_10 <- data.frame()

for (i in path){
  go_path_i <- go_path_id[go_path_id$go_path == i,]
  if (nrow(go_path_i) > 2){
    go_path_i <- go_path_i[1:2,]
  }
  go_path_10 <- rbind(go_path_10,go_path_i)
}
go_id <- go_path_10

counts <- read.csv("counts_cv.csv")
go_path_counts <- counts %>% filter(circRNA_ID %in% go_id$circRNA_ID)
write.csv(go_path_counts,"go_path_counts.csv",row.names = FALSE)

go_old <- read.csv('go_path_counts_train+test.csv')
go_corr <- read.csv('go_path_spearman.csv')
go_corr <- go_corr[order(go_corr$gene_correlations, decreasing = TRUE), ]
go_new <- go_corr[1:51,]
go_new <- go_new[order(go_new$go_path, decreasing = FALSE), ]
write.csv(go_new, "go_path_51.csv",row.names = FALSE)
go_new_count <- go_old %>% filter(circRNA_ID %in% go_new$circRNA_ID)
write.csv(go_new_count, "go_path_counts_train+test.csv",row.names = FALSE)

go <- read.csv('go_path_51.csv')
counts <- read.csv('counts_cv.csv')
go_new_count_cv <- counts %>% filter(circRNA_ID %in% go$circRNA_ID)
write.csv(go_new_count_cv, "go_path_counts.csv",row.names = FALSE)


#### Panel-3DG-cluster ####
go_path_id <- read.csv("go_path_id.csv")
info <- read_xlsx("3D_info.xlsx")
path <- unique(info$go_path)
k<- info$k
j <- 1
genome_path <- data.frame()

for (i in path){
  genome_path_i <- go_path_id[go_path_id$go_path == i,]
  num <- k[j]
  if (nrow(genome_path_i) > num){
    genome_path_i <- genome_path_i[1:num,]
  }
  genome_path <- rbind(genome_path,genome_path_i)
  j <- j+1
}

counts <- read.csv("counts_cv.csv")
genome_path_counts <- counts %>% filter(circRNA_ID %in% genome_path$circRNA_ID)
write.csv(genome_path_counts,"circ_3d_counts.csv",row.names = FALSE)


#### Panel-3DG-Radius ####
dis <- read.csv("all_bin_co_with_distance.csv")
bin_info <- read_xlsx("bin_info.xlsx")

bin_sign <- data.frame()
for (j in 0:17){
  start <- unlist(bin_info[(j+1), "Start_r"])[1]
  end <- unlist(bin_info[(j+1), "End_r"])[1]
  
  bin_1 <- data.frame()
  for (i in 1:nrow(dis)) {
    value <- unlist(dis[i, "distance"])[1]
    if (j == 17){
      if (value > start) {
        bin_l <- dis[i,]
        bin_1 <- rbind(bin_1, bin_l)
      }
    }
    else {
      if (value > start & value <= end) {
        bin_l <- dis[i,]
        bin_1 <- rbind(bin_1, bin_l)
      }
    }
  }
  bin_1$cobin <- j
  print(nrow(bin_1))
  bin_sign <- rbind(bin_sign, bin_1)
}
write.csv(bin_sign,"bin_sign.csv",row.names = FALSE)


## Panel-3DG-Radius1
bin_sign <- read.csv("bin_sign.csv")
bin_1 <- bin_sign[bin_sign$cobin == 1,]
radius_counts <- counts %>% filter(circRNA_ID %in% bin_1$circRNA_ID)
write.csv(radius_counts,"circ_radius_1_counts.csv",row.names = FALSE)


## Panel-3DG-Radius5
corr <- read.csv("circRNA_spearman_abs.csv")
bin_sign <- read.csv("bin_sign.csv")
bin_5 <- bin_sign[bin_sign$cobin == 5,]

corr_5 <- corr %>% filter(circRNA_ID %in% bin_5$circRNA_ID)
corr_5 <- corr_5[order(corr_5$gene_correlations, decreasing = TRUE), ]
circ_r <- corr_5[1:N,]
radius_counts <- counts %>% filter(circRNA_ID %in% circ_r$circRNA_ID)
write.csv(radius_counts,"circ_radius_5_counts.csv",row.names = FALSE)