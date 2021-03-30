require(vegan)

# Read files #
asvtab3 <- read_delim(here("data", "asvtab3_merged_subsamp_readnum.txt"), delim = "\t")
asvtab3c <- read_delim(here("data", "asvtab3c_merged_subsamp_pa.txt"), delim = "\t")
asvtab3_num <- asvtab3 %>% purrr::keep(is.numeric)
asvtab3c_num <- asvtab3c %>% purrr::keep(is.numeric)
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

#### Calculate slopes of the rarefaction curves at the endpoint ####
rarslopes <- c(NULL)
for (i in 1:dim(asvtab3_num)[2]) {
  rarslopes[i] <- rareslope(t(asvtab3_num[,i]), sample=sum(asvtab3_num[,i])-1)  
}
names(rarslopes) <- names(asvtab3_num)

#### Check correlation between ASV richness and slope of rarefaction curve ####
nasvs_sample_sizefract <- tibble("sample_sizefract" = names(asvtab3c_num), "nasvs" = colSums(asvtab3c_num))
slope_sample_sizefract <- tibble("sample_sizefract" = names(rarslopes), "slope" = rarslopes)

nasvs_slope <- left_join(nasvs_sample_sizefract, slope_sample_sizefract, by = "sample_sizefract")

cor.test(nasvs_slope$nasvs, nasvs_slope$slope) # cor = -0.12, p.value = 0.11 -> no correlation

#### Create rarefaction curves #### 
rareobj <- rarecurve(t(asvtab3_num),step=1000) 

str(rareobj)

#### Create new object for plotting in ggplot ####
rareo_msamp <- list()
for (i in 1:length(rareobj)) {
    rareo_msamp[[i]] <- cbind.data.frame("nASV" = rareobj[[i]], "step" = attr(rareobj[[i]], "Subsample"), "sample_sizefract" =rep(names(asvtab3_num)[i], length(rareobj[[i]])))
  }
  
rareo_msamp_tab <- bind_rows(rareo_msamp, .id = "column_label")
  
rareo_msamp_tab_meta <- left_join(rareo_msamp_tab, meta_tab_sample_sf, by = "sample_sizefract")
rareo_msamp_tab_meta <- rareo_msamp_tab_meta %>% mutate_at(vars(size_fraction), list( ~factor(., levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T)))
  
ggplot(rareo_msamp_tab_meta, aes(x = step, y = nASV, group = sample_sizefract))+
  geom_line(aes(color=size_fraction))+
  xlab("Number of reads")+
  ylab("Number of ASVs")
  