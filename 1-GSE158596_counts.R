###### Generate Counts Matrix Based on CIRI2 Results for GSE158596 ######

#---------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(conflicted)
library(readxl)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')

dir_path <- "GSE158596"
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

results <- matrix(0, nrow = length(all_ID), ncol = 94)
for (i in 1:94) {
  counts <- sapply(all_ID, function(x){sum(circRNA_ID[[i]] == x)})
  results[, i] <- counts
}

results_2 <- as_tibble(cbind(results, rowSums(results), index = 1:length(all_ID))) %>%
  filter(V95 >= 3) %>%
  dplyr::select(index) 

all_ID_index <- as_tibble(cbind(all_ID, index = 1:length(all_ID)))
all_ID_index$index <- as.double(all_ID_index$index)
selected_ID <- inner_join(all_ID_index, results_2, by = 'index')
save(selected_ID, file = 'Selected_ID.RData')
write.csv(selected_ID, '3Sample_circRNA_ID.csv', row.names = FALSE)


#----------
library(tidyverse)
library(ggplot2)
library(conflicted)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')
load(file = 'Selected_ID.RData')

##### gene_id
# inv16 -13
cancer_1 <- read_excel('ciri760.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_2 <- read_excel('ciri761.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_3 <- read_excel('ciri762.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_4 <- read_excel('ciri763.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_5 <- read_excel('ciri764.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_6 <- read_excel('ciri765.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_7 <- read_excel('ciri766.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_8 <- read_excel('ciri767.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_9 <- read_excel('ciri768.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_10 <- read_excel('ciri769.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_11 <- read_excel('ciri770.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_12 <- read_excel('ciri771.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_13 <- read_excel('ciri772.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# t821 -17
cancer_14 <- read_excel('ciri773.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_15 <- read_excel('ciri774.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_16 <- read_excel('ciri775.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_17 <- read_excel('ciri776.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_18 <- read_excel('ciri777.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_19 <- read_excel('ciri778.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_20 <- read_excel('ciri780.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_21 <- read_excel('ciri782.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_22 <- read_excel('ciri784.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_23 <- read_excel('ciri786.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_24 <- read_excel('ciri788.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_25 <- read_excel('ciri789.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_26 <- read_excel('ciri791.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_27 <- read_excel('ciri792.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_28 <- read_excel('ciri793.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_29 <- read_excel('ciri794.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_30 <- read_excel('ciri795.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_31 <- read_excel('ciri796.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# NPM1 -29+2
cancer_32 <- read_excel('ciri697.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_33 <- read_excel('ciri698.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_34 <- read_excel('ciri699.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_35 <- read_excel('ciri700.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_36 <- read_excel('ciri701.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_37 <- read_excel('ciri702.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_38 <- read_excel('ciri703.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_39 <- read_excel('ciri704.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_40 <- read_excel('ciri705.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_41 <- read_excel('ciri706.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_42 <- read_excel('ciri707.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_43 <- read_excel('ciri708.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_44 <- read_excel('ciri709.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_45 <- read_excel('ciri710.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_46 <- read_excel('ciri711.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_47 <- read_excel('ciri712.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_48 <- read_excel('ciri713.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_49 <- read_excel('ciri714.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_50 <- read_excel('ciri715.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_51 <- read_excel('ciri716.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_52 <- read_excel('ciri717.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_53 <- read_excel('ciri718.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_54 <- read_excel('ciri719.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_55 <- read_excel('ciri720.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_56 <- read_excel('ciri721.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_57 <- read_excel('ciri722.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_58 <- read_excel('ciri723.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_59 <- read_excel('ciri724.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_60 <- read_excel('ciri725.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_61 <- read_excel('ciri726.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_62 <- read_excel('ciri727.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# SF3B1 -3+1
cancer_63 <- read_excel('ciri728.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_64 <- read_excel('ciri729.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_65 <- read_excel('ciri730.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_66 <- read_excel('ciri731.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# SRFR2 -5
cancer_67 <- read_excel('ciri732.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_68 <- read_excel('ciri733.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_69 <- read_excel('ciri734.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_70 <- read_excel('ciri735.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_71 <- read_excel('ciri736.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# SF382_SRFR2 -+1
cancer_72 <- read_excel('ciri737.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# U2AF1 -6
cancer_73 <- read_excel('ciri738.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_74 <- read_excel('ciri739.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_75 <- read_excel('ciri740.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_76 <- read_excel('ciri741.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_77 <- read_excel('ciri742.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
cancer_78 <- read_excel('ciri743.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

# non-cancer
notcancer_1 <- read_excel('ciri754.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_2 <- read_excel('ciri755.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_3 <- read_excel('ciri756.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_4 <- read_excel('ciri757.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_5 <- read_excel('ciri758.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_6 <- read_excel('ciri759.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)

notcancer_7 <- read_excel('ciri744.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_8 <- read_excel('ciri745.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_9 <- read_excel('ciri746.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_10 <- read_excel('ciri747.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_11 <- read_excel('ciri748.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_12 <- read_excel('ciri749.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_13 <- read_excel('ciri750.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_14 <- read_excel('ciri751.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_15 <- read_excel('ciri752.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)
notcancer_16 <- read_excel('ciri753.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID)


notcancer <- bind_rows(notcancer_1, notcancer_2, notcancer_3, notcancer_4, notcancer_5, notcancer_6, notcancer_7, notcancer_8, 
                       notcancer_9, notcancer_10, notcancer_11, notcancer_12, notcancer_13, notcancer_14, 
                       notcancer_15, notcancer_16,
                       .id = 'notcancer')
cancer <- bind_rows(cancer_1, cancer_2, cancer_3, cancer_4, cancer_5, cancer_6, 
                    cancer_7, cancer_8, cancer_9, cancer_10,cancer_11, cancer_12, cancer_13, cancer_14, cancer_15, cancer_16, 
                    cancer_17, cancer_18, cancer_19, cancer_20,cancer_21, cancer_22, cancer_23, cancer_24, cancer_25, cancer_26, 
                    cancer_27, cancer_28, cancer_29, cancer_30,cancer_31, cancer_32, cancer_33, cancer_34, cancer_35, cancer_36, 
                    cancer_37, cancer_38, cancer_39, cancer_40,cancer_41, cancer_42, cancer_43, cancer_44, cancer_45, cancer_46, 
                    cancer_47, cancer_48, cancer_49, cancer_50,cancer_51, cancer_52, cancer_53, cancer_54, cancer_55, cancer_56, 
                    cancer_57, cancer_58, cancer_59, cancer_60,cancer_61, cancer_62, cancer_63, cancer_64, cancer_65, cancer_66, 
                    cancer_67, cancer_68, cancer_69, cancer_70,cancer_71, cancer_72, cancer_73, cancer_74, cancer_75, cancer_76, 
                    cancer_77, cancer_78,
                    .id = 'cancer')
data <- bind_rows(notcancer, cancer)
write.csv(data, 'GSE158596_data.csv', row.names = FALSE)

gene_id <- bind_rows(notcancer, cancer) %>%
  dplyr::select(circRNA_ID, gene_id) %>%
  dplyr::filter(gene_id != 'n/a')

max_comma_count <- max(str_count(gene_id$gene_id, ","))
cat("The maximum comma count:", max_comma_count, "\n")
library(tidyr)
gene_id <- separate(gene_id, gene_id, into = paste0("gene_", 116), sep = ",")
gene_id_long <- pivot_longer(gene_id, cols = starts_with("gene_"), names_to = "gene_column", values_to = "gene_id") %>%
  dplyr::filter(!is.na(gene_id)) %>%
  dplyr::filter(gene_id != '') %>%
  select(circRNA_ID, gene_id) %>%
  distinct()

write.csv(gene_id_long, 'GSE158596_gene_id.csv', row.names = FALSE)

##### merge circRNA data #####
library(tidyverse)
library(ggplot2)
library(conflicted)
conflict_prefer('filter', 'dplyr')
conflict_prefer('select', 'dplyr')
load(file = 'Selected_ID.RData')

# inv16 -13
cancer_1 <- read_excel('ciri760.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_2 <- read_excel('ciri761.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_3 <- read_excel('ciri762.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_4 <- read_excel('ciri763.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_5 <- read_excel('ciri764.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_6 <- read_excel('ciri765.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_7 <- read_excel('ciri766.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_8 <- read_excel('ciri767.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_9 <- read_excel('ciri768.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_10 <- read_excel('ciri769.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_11 <- read_excel('ciri770.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_12 <- read_excel('ciri771.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_13 <- read_excel('ciri772.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# t821 -17
cancer_14 <- read_excel('ciri773.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_15 <- read_excel('ciri774.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_16 <- read_excel('ciri775.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_17 <- read_excel('ciri776.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_18 <- read_excel('ciri777.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_19 <- read_excel('ciri778.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_20 <- read_excel('ciri780.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_21 <- read_excel('ciri782.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_22 <- read_excel('ciri784.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_23 <- read_excel('ciri786.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_24 <- read_excel('ciri788.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_25 <- read_excel('ciri789.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_26 <- read_excel('ciri791.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_27 <- read_excel('ciri792.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_28 <- read_excel('ciri793.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_29 <- read_excel('ciri794.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_30 <- read_excel('ciri795.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_31 <- read_excel('ciri796.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# NPM1 -29+2
cancer_32 <- read_excel('ciri697.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_33 <- read_excel('ciri698.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_34 <- read_excel('ciri699.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_35 <- read_excel('ciri700.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_36 <- read_excel('ciri701.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_37 <- read_excel('ciri702.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_38 <- read_excel('ciri703.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_39 <- read_excel('ciri704.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_40 <- read_excel('ciri705.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_41 <- read_excel('ciri706.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_42 <- read_excel('ciri707.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_43 <- read_excel('ciri708.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_44 <- read_excel('ciri709.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_45 <- read_excel('ciri710.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_46 <- read_excel('ciri711.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_47 <- read_excel('ciri712.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_48 <- read_excel('ciri713.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_49 <- read_excel('ciri714.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_50 <- read_excel('ciri715.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_51 <- read_excel('ciri716.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_52 <- read_excel('ciri717.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_53 <- read_excel('ciri718.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_54 <- read_excel('ciri719.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_55 <- read_excel('ciri720.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_56 <- read_excel('ciri721.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_57 <- read_excel('ciri722.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_58 <- read_excel('ciri723.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_59 <- read_excel('ciri724.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_60 <- read_excel('ciri725.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_61 <- read_excel('ciri726.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_62 <- read_excel('ciri727.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# SF3B1 -3+1
cancer_63 <- read_excel('ciri728.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_64 <- read_excel('ciri729.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_65 <- read_excel('ciri730.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_66 <- read_excel('ciri731.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# SRFR2 -5
cancer_67 <- read_excel('ciri732.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_68 <- read_excel('ciri733.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_69 <- read_excel('ciri734.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_70 <- read_excel('ciri735.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_71 <- read_excel('ciri736.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# SF382_SRFR2 -+1
cancer_72 <- read_excel('ciri737.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# U2AF1 -6
cancer_73 <- read_excel('ciri738.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_74 <- read_excel('ciri739.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_75 <- read_excel('ciri740.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_76 <- read_excel('ciri741.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_77 <- read_excel('ciri742.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
cancer_78 <- read_excel('ciri743.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

# non-cancer
notcancer_1 <- read_excel('ciri754.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_2 <- read_excel('ciri755.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_3 <- read_excel('ciri756.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_4 <- read_excel('ciri757.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_5 <- read_excel('ciri758.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_6 <- read_excel('ciri759.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

notcancer_7 <- read_excel('ciri744.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_8 <- read_excel('ciri745.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_9 <- read_excel('ciri746.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_10 <- read_excel('ciri747.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_11 <- read_excel('ciri748.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_12 <- read_excel('ciri749.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_13 <- read_excel('ciri750.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_14 <- read_excel('ciri751.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_15 <- read_excel('ciri752.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)
notcancer_16 <- read_excel('ciri753.xlsx') %>%
  filter(circRNA_ID %in% selected_ID$all_ID) %>%
  select(circRNA_ID, `#junction_reads`)

notcancer <- bind_rows(notcancer_1, notcancer_2, notcancer_3, notcancer_4, notcancer_5, notcancer_6, notcancer_7, notcancer_8, 
                       notcancer_9, notcancer_10, notcancer_11, notcancer_12, notcancer_13, notcancer_14, 
                       notcancer_15, notcancer_16,
                       .id = 'notcancer') %>%
  pivot_wider(names_from = notcancer, values_from = `#junction_reads`)
notcancer[is.na(notcancer)] <- 0
colnames(notcancer) <- c('circRNA_ID', paste0('notcancer', '_', 1:16))

cancer <- bind_rows(cancer_1, cancer_2, cancer_3, cancer_4, cancer_5, cancer_6, 
                    cancer_7, cancer_8, cancer_9, cancer_10,cancer_11, cancer_12, cancer_13, cancer_14, cancer_15, cancer_16, 
                    cancer_17, cancer_18, cancer_19, cancer_20,cancer_21, cancer_22, cancer_23, cancer_24, cancer_25, cancer_26, 
                    cancer_27, cancer_28, cancer_29, cancer_30,cancer_31, cancer_32, cancer_33, cancer_34, cancer_35, cancer_36, 
                    cancer_37, cancer_38, cancer_39, cancer_40,cancer_41, cancer_42, cancer_43, cancer_44, cancer_45, cancer_46, 
                    cancer_47, cancer_48, cancer_49, cancer_50,cancer_51, cancer_52, cancer_53, cancer_54, cancer_55, cancer_56, 
                    cancer_57, cancer_58, cancer_59, cancer_60,cancer_61, cancer_62, cancer_63, cancer_64, cancer_65, cancer_66, 
                    cancer_67, cancer_68, cancer_69, cancer_70,cancer_71, cancer_72, cancer_73, cancer_74, cancer_75, cancer_76, 
                    cancer_77, cancer_78,
                    .id = 'cancer') %>%
  pivot_wider(names_from = cancer, values_from = `#junction_reads`)
cancer[is.na(cancer)] <- 0
colnames(cancer) <- c('circRNA_ID', paste0('cancer', '_', 1:78))

counts <- full_join(notcancer, cancer, by = 'circRNA_ID')
counts[is.na(counts)] <- 0
counts <- as.data.frame(counts)
write.csv(counts, "GSE158596_counts.csv", row.names = FALSE)