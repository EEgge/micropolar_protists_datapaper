library(here)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(purrr)
library(dplyr)
library(readr)
library(readxl)


#### Read files ####
asvtab_wtax <- read_excel(here("metapr2/export", "metapr2_wide_asv_set_207_208_209_Eukaryota_nodups.xlsx"))
dim(asvtab_wtax)
length(unique(asvtab_wtax$asv_code))
#7054 x 208
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]
env_tab <- read_delim(here("data", "meta_cleanMP.txt"), delim = "\t")


asvtab_wtax <- asvtab_wtax %>% add_column(divisionlong = NA)
for (i in 1:dim(asvtab_wtax)[1]) {
  if (asvtab_wtax$division[i] %in% c("Dinoflagellata", "Ochrophyta")) {
    asvtab_wtax$divisionlong[i] <- paste(asvtab_wtax$division[i], asvtab_wtax$class[i], sep = "_")
  } else {
    asvtab_wtax$divisionlong[i] <- asvtab_wtax$division[i]
  }
}

asvtab_wo_metazoa_embr <- asvtab_wtax %>% filter(division !="Metazoa", class != "Embryophyceae")
dim(asvtab_wo_metazoa_embr) #6545 x 209
length(unique(asvtab_wo_metazoa_embr$asv_code))
write_delim(asvtab_wo_metazoa_embr, "data/asvtab1_nonmerged_readnum.txt", delim = "\t")



asvtab_num <- asvtab_wo_metazoa_embr %>% select_if(is.numeric)
asvtab_notnum <- asvtab_wo_metazoa_embr %>% select_if(negate(is.numeric))

asvtab_pa <- asvtab_num
asvtab_pa[asvtab_pa>0] <- 1
asvtab_nonsubsamp_pa <- bind_cols(asvtab_notnum, asvtab_pa)
write_delim(asvtab_nonsubsamp_pa, "data/asvtab3_nonmerged_pa.txt", delim = "\t")


colSums(asvtab_num)

meta_tab_nodups <- meta_tab[!duplicated(meta_tab$sample_sizefract),]

asvtab_wo_metazoa_embr <- read_delim(here("data", "asvtab1_nonmerged_readnum.txt"), delim = "\t")
asvtab_num <- asvtab_wo_metazoa_embr %>% purrr::keep(is.numeric)
#### pivot longer, then merge "sample_sizefract" with replicates:  ####
asvtab_num_pivot <- pivot_longer(asvtab_wo_metazoa_embr, cols=colnames(asvtab_num), names_to = "seq_event")
asvtab_num_pivot_sampsf <- left_join(asvtab_num_pivot, meta_tab_seq_event, by = "seq_event")
asvtab_num_pivot_sampsf_merge <- asvtab_num_pivot_sampsf %>% group_by(sample_sizefract, asv_code, kingdom, supergroup, division, class, family, order, genus, species, divisionlong, sequence) %>% 
                                 select(-file_name) %>% summarise_if(is.numeric, sum)

#pivot wider prior to rrarefy:
merged_asv_tab <- pivot_wider(asvtab_num_pivot_sampsf_merge, id_cols = c(asv_code, kingdom, supergroup, division, divisionlong, class, family, order, genus, species, sequence), names_from = sample_sizefract, values_from = value)
dim(merged_asv_tab) 
write_delim(merged_asv_tab, "data/asvtab1b_merged_readnum.txt", delim = "\t")

#then subsample even read number per sample_sizefract
#number of reads in each sample_sizefract:
nreads <- merged_asv_tab %>% purrr::keep(is.numeric) %>% colSums()

colsumtab <- cbind.data.frame("sample_sizefract" = names(nreads), "totreads" = nreads)
left_join(colsumtab, meta_tab_nodups, by = "sample_sizefract") %>% filter(fraction_max == 200) %>% 
  select(sample_sizefract, totreads)


#### Subsample to equal read number within size fractions ####
#Subsample each size fraction separately:

perm_rr <- 100
merged_asv_0.4_3 <- merged_asv_tab %>% ungroup() %>% select(meta_tab %>% filter(fraction_min == 0.4, fraction_max == 3) %>% pull(sample_sizefract))
merged_asv_0.4_3_df <- as.data.frame(merged_asv_0.4_3, row.names = merged_asv_tab$asv_code)
dim(merged_asv_0.4_3_df) #6536 x 44

rrarearray_0.4_3_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_0.4_3_df)), ncol(t(merged_asv_0.4_3_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_0.4_3_40000[,,i] <- rrarefy(t(merged_asv_0.4_3_df), 40000)
}

asvtab_rravg_0.4_3_40000 <- as.data.frame(round(t(apply(rrarearray_0.4_3_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_0.4_3_40000)
colnames(asvtab_rravg_0.4_3_40000) <- names(merged_asv_0.4_3_df)
asvtab_rravg_0.4_3_40000$asv_code <- merged_asv_tab$asv_code
dim(asvtab_rravg_0.4_3_40000)

#### rarefy 3_10 ####
merged_asv_3_10 <- merged_asv_tab %>% ungroup() %>% select(meta_tab %>% filter(fraction_min == 3, fraction_max == 10) %>% pull(sample_sizefract))
merged_asv_3_10_df <- as.data.frame(merged_asv_3_10, row.names = merged_asv_tab$asv_code)
dim(merged_asv_3_10_df)

rrarearray_3_10_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_3_10_df)), ncol(t(merged_asv_3_10_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_3_10_40000[,,i] <- rrarefy(t(merged_asv_3_10_df), 40000)
}

asvtab_rravg_3_10_40000 <- as.data.frame(round(t(apply(rrarearray_3_10_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_3_10_40000)
colnames(asvtab_rravg_3_10_40000) <- names(merged_asv_3_10_df)
asvtab_rravg_3_10_40000$asv_code <- merged_asv_tab$asv_code
dim(asvtab_rravg_3_10_40000)

#### Rarefy 10_50 ####
merged_asv_10_50 <- merged_asv_tab %>% ungroup() %>% select(meta_tab %>% filter(fraction_min == 10, fraction_max == 50) %>% pull(sample_sizefract))
merged_asv_10_50_df <- as.data.frame(merged_asv_10_50, row.names = merged_asv_tab$asv_code)
dim(merged_asv_10_50_df)

rrarearray_10_50_40000 <- array(data=NA, dim=c(nrow(t(merged_asv_10_50_df)), ncol(t(merged_asv_10_50_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_10_50_40000[,,i] <- rrarefy(t(merged_asv_10_50_df), 40000)
}

asvtab_rravg_10_50_40000 <- as.data.frame(round(t(apply(rrarearray_10_50_40000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_10_50_40000)
colnames(asvtab_rravg_10_50_40000) <- names(merged_asv_10_50_df)
asvtab_rravg_10_50_40000$asv_code <- merged_asv_tab$asv_code
dim(asvtab_rravg_10_50_40000)

#### rarefy 50_200 ####
merged_asv_50_200 <- merged_asv_tab %>% ungroup() %>% select(meta_tab %>% filter(fraction_min == 50, fraction_max == 200) %>% pull(sample_sizefract))
merged_asv_50_200_df <- as.data.frame(merged_asv_50_200, row.names = merged_asv_tab$asv_code)
dim(merged_asv_50_200_df)

rrarearray_50_200_8000 <- array(data=NA, dim=c(nrow(t(merged_asv_50_200_df)), ncol(t(merged_asv_50_200_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_50_200_8000[,,i] <- rrarefy(t(merged_asv_50_200_df), 8000)
}

asvtab_rravg_50_200_8000 <- as.data.frame(round(t(apply(rrarearray_50_200_8000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_50_200_8000)
colnames(asvtab_rravg_50_200_8000) <- names(merged_asv_50_200_df)
asvtab_rravg_50_200_8000$asv_code <- merged_asv_tab$asv_code
dim(asvtab_rravg_50_200_8000)

#### rrarefy 3-180 ####
merged_asv_3_180 <- merged_asv_tab %>% ungroup() %>% select(meta_tab %>% filter(fraction_min == 3, fraction_max == 180) %>% pull(sample_sizefract))
merged_asv_3_180_df <- as.data.frame(merged_asv_3_180, row.names = merged_asv_tab$asv_code)

rrarearray_3_180_88000 <- array(data=NA, dim=c(nrow(t(merged_asv_3_180_df)), ncol(t(merged_asv_3_180_df)), perm_rr))
for (i in 1:perm_rr) {
  rrarearray_3_180_88000[,,i] <- rrarefy(t(merged_asv_3_180_df), 88000)
}

asvtab_rravg_3_180_88000 <- as.data.frame(round(t(apply(rrarearray_3_180_88000, c(1,2), mean)), digits = 0))

dim(asvtab_rravg_3_180_88000)
colnames(asvtab_rravg_3_180_88000) <- names(merged_asv_3_180_df)
asvtab_rravg_3_180_88000$asv_code <- merged_asv_tab$asv_code
dim(asvtab_rravg_3_180_88000)

#### join subsampled data frames ####
asv_tab_merged_subsamp <- asvtab_notnum %>% left_join(asvtab_rravg_0.4_3_40000, by = "asv_code") %>% 
  left_join(asvtab_rravg_3_10_40000, by = "asv_code") %>% 
  left_join(asvtab_rravg_10_50_40000, by = "asv_code") %>%
  left_join(asvtab_rravg_50_200_8000, by = "asv_code") %>%
  left_join(asvtab_rravg_3_180_88000, by = "asv_code")
dim(asv_tab_merged_subsamp)
write_delim(asv_tab_merged_subsamp, "data/asv_tab_merged_subsamp.txt", delim = "\t")

#### Create proportions after subsampling ####
asvtab_subsamp_num <- asv_tab_merged_subsamp %>% select_if(is.numeric)
asvtab_subsamp_notnum <- asv_tab_merged_subsamp %>% select_if(negate(is.numeric))
asvtab_subsamp_num_prop <- sweep(asvtab_subsamp_num, 2 , colSums(asvtab_subsamp_num), FUN = "/")

colSums(asvtab_subsamp_num_prop)
asvtab_subsamp_prop_tax <- bind_cols(asvtab_subsamp_notnum, asvtab_subsamp_num_prop)
write_delim(asvtab_subsamp_prop_tax, "data/asvtab_subsamp_prop_wtax.txt", delim = "\t")

dim(asvtab_subsamp_prop_tax) 


#### Create presence-absence after subsampling asv table ####
asvtab_subsamp_notnum <- asvtab_subsamp_prop_wtax %>% select_if(negate(is.numeric))
asvtab_subsamp_num_prop <- asvtab_subsamp_prop_wtax %>% select_if(is.numeric)
asvtab_subsamp_num_prop_pa <- asvtab_subsamp_num_prop
asvtab_subsamp_num_prop_pa[asvtab_subsamp_num_prop_pa >0] <- 1
asvtab_subsamp_pa_wtax <- bind_cols(asvtab_subsamp_notnum, asvtab_subsamp_num_prop_pa)
write_delim(asvtab_subsamp_pa_wtax, "data/asvtab_subsamp_pa_wtax.txt", delim = "\t")

