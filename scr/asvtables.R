library(here)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(purrr)
library(dplyr)
library(readr)
library(readxl)


# Read files #
asvtab_wtax <- read_excel(here("metapr2/export", "metapr2_wide_asv_set_207_208_209_Eukaryota.xlsx"))
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

# Add extra taxonomic levels to split Dinoflagellata and Ocrophyta into sub-groups
asvtab_wtax <- asvtab_wtax %>% add_column(divisionlong = NA)
for (i in 1:dim(asvtab_wtax)[1]) {
  if (asvtab_wtax$division[i] %in% c("Dinoflagellata", "Ochrophyta")) {
    asvtab_wtax$divisionlong[i] <- paste(asvtab_wtax$division[i], asvtab_wtax$class[i], sep = "_")
  } else {
    asvtab_wtax$divisionlong[i] <- asvtab_wtax$division[i]
  }
}

#### asvtab1 - Filter out Metazoa and Embryophyceae ####
asvtab_wo_metazoa_embr <- asvtab_wtax %>% filter(division !="Metazoa", class != "Embryophyceae")
dim(asvtab_wo_metazoa_embr)
write_delim(asvtab_wo_metazoa_embr, "data/asvtab1_nonmerged_readnum.txt", delim = "\t")

# Transform asvtab1 non-merged, non-subsampled asv table to proportions #
asvtab_num <- asvtab_wo_metazoa_embr %>% purrr::keep(is.numeric)
asvtab_notnum <- asvtab_wo_metazoa_embr %>% purrr::keep(negate(is.numeric))
asvtab_num_prop <- sweep(asvtab_num, 2 , colSums(asvtab_num), FUN = "/")

colSums(asvtab_num)
asvtab_prop_tax <- cbind.data.frame(asvtab_notnum, asvtab_num_prop)
write_delim(asvtab_prop_tax, "data/asvtab1b_nonmerged_prop.txt", delim = "\t")

# Transform asvtab1 non-merged, non-subsampled asv table to presence-absence #
asvtab_pa <- asvtab_num
asvtab_pa[asvtab_pa>0] <- 1
asvtab_nonsubsamp_pa <- bind_cols(asvtab_notnum, asvtab_pa)
write_delim(asvtab_nonsubsamp_pa, "data/asvtab1c_nonmerged_pa.txt", delim = "\t")


#### asvtab2 - Merge sequencing_events that represent the same sample_sizefract ####
asvtab_num_pivot <- pivot_longer(asvtab_wo_metazoa_embr, cols=colnames(asvtab_num), names_to = "seq_event")
asvtab_num_pivot_sampsf <- left_join(asvtab_num_pivot, meta_tab_seq_event, by = "seq_event")
asvtab_num_pivot_sampsf_merge <- asvtab_num_pivot_sampsf %>% group_by(sample_sizefract, asv_code, kingdom, supergroup, division, class, family, order, genus, species, divisionlong, sequence) %>% 
                                 select(-file_name) %>% summarise_if(is.numeric, sum)

asvtab_merged <- pivot_wider(asvtab_num_pivot_sampsf_merge, id_cols = c(asv_code, kingdom, supergroup, division, divisionlong, class, family, order, genus, species, sequence), names_from = sample_sizefract, values_from = value)

write_delim(asvtab_merged, "data/asvtab2_merged_readnum.txt", delim = "\t")

# Transform asvtab2 merged, non-subsampled asv table to proportions #
asvtab_merged_num <- asvtab_merged %>% purrr::keep(is.numeric)
asvtab_merged_notnum <- asvtab_merged %>% purrr::keep(negate(is.numeric))
asvtab_merged_prop <- sweep(asvtab_merged_num, 2 , colSums(asvtab_merged_num), FUN = "/")

asvtab_merged_prop_tax <- cbind.data.frame(asvtab_merged_notnum, asvtab_merged_prop)
write_delim(asvtab_prop_tax, "data/asvtab2b_merged_prop.txt", delim = "\t")

# Transform asvtab2 non-merged, non-subsampled asv table to presence-absence #
asvtab_merged_pa <- asvtab_merged_num
asvtab_merged_pa[asvtab_merged_pa>0] <- 1
asvtab_merged_pa_wtax <- bind_cols(asvtab_merged_notnum, asvtab_merged_pa)
write_delim(asvtab_merged_pa_wtax, "data/asvtab2c_merged_pa.txt", delim = "\t")


#### asvtab3 Subsample to equal read number within size fractions ####
#Subsample each size fraction separately:

perm_rr <- 100
# 0.45-3 #
merged_asv_0.4_3 <- asvtab_merged %>% ungroup() %>% select(meta_tab_sample_sf %>% filter(fraction_min == 0.4, fraction_max == 3) %>% pull(sample_sizefract))
merged_asv_0.4_3_df <- as.data.frame(merged_asv_0.4_3, row.names = asvtab_merged$asv_code)

rrarearray_0.4_3_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_0.4_3_df)), ncol(t(merged_asv_0.4_3_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_0.4_3_40000[,,i] <- rrarefy(t(merged_asv_0.4_3_df), 40000)
}

asvtab_rravg_0.4_3_40000 <- as.data.frame(round(t(apply(rrarearray_0.4_3_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_0.4_3_40000)
colnames(asvtab_rravg_0.4_3_40000) <- names(merged_asv_0.4_3_df)
asvtab_rravg_0.4_3_40000$asv_code <- asvtab_merged$asv_code
dim(asvtab_rravg_0.4_3_40000)

# 3-10 #
merged_asv_3_10 <- asvtab_merged %>% ungroup() %>% select(meta_tab_sample_sf %>% filter(fraction_min == 3, fraction_max == 10) %>% pull(sample_sizefract))
merged_asv_3_10_df <- as.data.frame(merged_asv_3_10, row.names = asvtab_merged$asv_code)
dim(merged_asv_3_10_df)

rrarearray_3_10_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_3_10_df)), ncol(t(merged_asv_3_10_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_3_10_40000[,,i] <- rrarefy(t(merged_asv_3_10_df), 40000)
}

asvtab_rravg_3_10_40000 <- as.data.frame(round(t(apply(rrarearray_3_10_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_3_10_40000)
colnames(asvtab_rravg_3_10_40000) <- names(merged_asv_3_10_df)
asvtab_rravg_3_10_40000$asv_code <- asvtab_merged$asv_code
dim(asvtab_rravg_3_10_40000)

# 10-50 #
merged_asv_10_50 <- asvtab_merged %>% ungroup() %>% select(meta_tab_sample_sf %>% filter(fraction_min == 10, fraction_max == 50) %>% pull(sample_sizefract))
merged_asv_10_50_df <- as.data.frame(merged_asv_10_50, row.names = asvtab_merged$asv_code)
dim(merged_asv_10_50_df)

rrarearray_10_50_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_10_50_df)), ncol(t(merged_asv_10_50_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_10_50_40000[,,i] <- rrarefy(t(merged_asv_10_50_df), 40000)
}

asvtab_rravg_10_50_40000 <- as.data.frame(round(t(apply(rrarearray_10_50_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_10_50_40000)
colnames(asvtab_rravg_10_50_40000) <- names(merged_asv_10_50_df)
asvtab_rravg_10_50_40000$asv_code <- asvtab_merged$asv_code
dim(asvtab_rravg_10_50_40000)

# 50-200 #
merged_asv_50_200 <- asvtab_merged %>% ungroup() %>% select(meta_tab_sample_sf %>% filter(fraction_min == 50, fraction_max == 200) %>% pull(sample_sizefract))
merged_asv_50_200_df <- as.data.frame(merged_asv_50_200, row.names = asvtab_merged$asv_code)
dim(merged_asv_50_200_df)

rrarearray_50_200_8000 <- array(data=NA, dim=c(nrow(t(merged_asv_50_200_df)), ncol(t(merged_asv_50_200_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_50_200_8000[,,i] <- rrarefy(t(merged_asv_50_200_df), 8000)
}

asvtab_rravg_50_200_8000 <- as.data.frame(round(t(apply(rrarearray_50_200_8000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_50_200_8000)
colnames(asvtab_rravg_50_200_8000) <- names(merged_asv_50_200_df)
asvtab_rravg_50_200_8000$asv_code <- asvtab_merged$asv_code
dim(asvtab_rravg_50_200_8000)

# 3-180 #
merged_asv_3_180 <- asvtab_merged %>% ungroup() %>% select(meta_tab_sample_sf %>% filter(fraction_min == 3, fraction_max == 180) %>% pull(sample_sizefract))
merged_asv_3_180_df <- as.data.frame(merged_asv_3_180, row.names = asvtab_merged$asv_code)

rrarearray_3_180_88000 <- array(data=NA, dim=c(nrow(t(merged_asv_3_180_df)), ncol(t(merged_asv_3_180_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_3_180_88000[,,i] <- rrarefy(t(merged_asv_3_180_df), 88000)
}

asvtab_rravg_3_180_88000 <- as.data.frame(round(t(apply(rrarearray_3_180_88000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_3_180_88000)
colnames(asvtab_rravg_3_180_88000) <- names(merged_asv_3_180_df)
asvtab_rravg_3_180_88000$asv_code <- asvtab_merged$asv_code
dim(asvtab_rravg_3_180_88000)

# join subsampled data frames #
asv_tab_merged_subsamp <- asvtab_notnum %>% left_join(asvtab_rravg_0.4_3_40000, by = "asv_code") %>% 
  left_join(asvtab_rravg_3_10_40000, by = "asv_code") %>% 
  left_join(asvtab_rravg_10_50_40000, by = "asv_code") %>%
  left_join(asvtab_rravg_50_200_8000, by = "asv_code") %>%
  left_join(asvtab_rravg_3_180_88000, by = "asv_code")
dim(asv_tab_merged_subsamp)
write_delim(asv_tab_merged_subsamp, "data/asvtab3_merged_subsamp_readnum.txt", delim = "\t")

# Transform asvtab3 merged, subsampled to proportions #
asvtab_subsamp_num <- asv_tab_merged_subsamp %>% purrr::keep(is.numeric)
asvtab_subsamp_notnum <- asv_tab_merged_subsamp %>% purrr::keep(negate(is.numeric))
asvtab_subsamp_num_prop <- sweep(asvtab_subsamp_num, 2 , colSums(asvtab_subsamp_num), FUN = "/")

colSums(asvtab_subsamp_num_prop)
asvtab_subsamp_prop_tax <- bind_cols(asvtab_subsamp_notnum, asvtab_subsamp_num_prop)
write_delim(asvtab_subsamp_prop_tax, "data/asvtab3b_merged_subsamp_prop.txt", delim = "\t")

# Transform asvtab3 merged, subsampled to presence-absence #
asvtab_subsamp_notnum <- asvtab_subsamp_prop_tax %>% purrr::keep(negate(is.numeric))
asvtab_subsamp_num_prop <- asvtab_subsamp_prop_tax %>% purrr::keep(is.numeric)
asvtab_subsamp_num_pa <- asvtab_subsamp_num_prop
asvtab_subsamp_num_pa[asvtab_subsamp_num_pa >0] <- 1
asvtab_subsamp_pa_wtax <- bind_cols(asvtab_subsamp_notnum, asvtab_subsamp_num_prop_pa)
write_delim(asvtab_subsamp_pa_wtax, "data/asvtab3c_merged_subsamp_pa.txt", delim = "\t")

