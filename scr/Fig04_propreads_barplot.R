library(shiny)
library(here)
library(ggplot2)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)
library(mvabund)
library(RColorBrewer)
library(reshape2)
library(purrr)
library(dplyr)
library(plotly)
library(indicspecies)
library(adespatial)
library(fossil)
library(ggdendro)
library(dendextend)
library(readr)
library(readxl)
library(xtable)

#### Read files ####
# Meta data # 
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

# asv tables #
taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")
asvtab5 <- read_delim(here("data", "asvtab5_merged_subsamp_prop.txt"), delim = "\t")


##### BARPLOT #####
#### Make legend ####
bpcol <- read_delim(here("data","Fig04_colors.txt"), col_names = T, delim = "\t", comment = "")
bpcol <- bpcol %>% mutate_at(vars(Divisionlong), list( ~factor(., levels = unique(bpcol$Divisionlong), ordered = T))) %>% filter(Divisionlong != "Centroheliozoa")
bpcolvec <- as.character(bpcol$col)
names(bpcolvec) <- bpcol$Divisionlong

bp_leg_plot <- ggplot(bpcol, aes(x = Divisionlong, fill = Divisionlong))+
  geom_bar()+
  scale_fill_manual(values = bpcolvec)+
  theme(legend.position = "right")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)),
         fill = guide_legend(nrow = 4, byrow = F))+
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.7, "lines"))+
  labs(fill = "Division/Class")
bp_leg_plot
bp_legnd <- get_legend(bp_leg_plot)




#### Group by Division ####
groupby_division <- asvtab5 %>% group_by(divisionlong) %>% summarise_if(is.numeric, sum)

#### Function to identify taxa that were present with < 5% of the reads in all samples ####
limfun <- function(x) {
  ifelse(x>=0.05,1,0)
}

## Transform grouped df to matrix, to be able to apply limfun ##
groupby_mat <- as.matrix(groupby_division[,-1, drop = FALSE])
groupby_bin <- apply(groupby_mat, 2, limfun)

if (is.vector(groupby_bin)) {
  groupby_bin2 <- as.data.frame(as.list(groupby_bin))
} else {
  groupby_bin2 <- as.data.frame(groupby_bin)
  }


#### Identify taxa present as < 5% in all samples (low prop. taxa) ####
taxgroupspre_bin_sums <- groupby_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
taxgroupspre_bin_sums$taxgroups <- groupby_division %>% pull(get("divisionlong"))
taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% select(taxgroups)

taxgroups_select <- groupby_division %>% filter(divisionlong %in% taxgroupspre_bin_sums_yes$taxgroups)

taxgroups_other <- groupby_division %>% filter(!divisionlong %in% taxgroupspre_bin_sums_yes$taxgroups)


## If low prop. taxa detected, sum read number in each sample ##
if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
  othersum <- colSums(taxgroups_other[,-1])
} else {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
  othersum <- NULL
}

#### Create data frame with >= 5% taxa, and all < 5% taxa grouped together ####
taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)

#### Pivot table ####
division_pivot <- pivot_longer(taxgroups_select3, cols = !Taxonomic_group, names_to = "sample_sizefract")

division_pivot <- division_pivot %>% mutate(Taxonomic_group = factor(Taxonomic_group, levels=.env$Taxonomic_group, ordered = T))


#### Join pivoted taxa table with metadata ####
division_pivot_wmeta <- left_join(division_pivot, meta_tab_sample_sf, by = "sample_sizefract")

division_pivot_wmeta_merge <- division_pivot_wmeta %>% group_by(Taxonomic_group, sample_sizefract, station_depth, env_sample, size_fraction, month, station, collection_method) %>%
                              summarise_if(is.numeric, mean) %>% 
                              select(-dna_concentration, -`260_280`, -`260_230`, -sample_sequence)

division_pivot_wmeta2 <- division_pivot_wmeta_merge %>% mutate(station_depth = factor(station_depth, levels = unique(meta_tab_sample_sf$station_depth), ordered = T),
                                                          month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T), 
                                                          collection_method = factor(collection_method, levels = c("niskin", "net")),
                                                          size_fraction = factor(size_fraction, levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T))

division_pivot_wmeta3 <- division_pivot_wmeta2 %>% mutate(size_fraction=recode(size_fraction, "0.4_3" = "0.45-3", "3_10" = "3-10", "10_50" = "10-50", "50_200" = "50-200", "3_180" = "3-180"))
division_pivot_wmeta4 <- division_pivot_wmeta3 %>% mutate(newsizefract = forcats::fct_collapse(size_fraction, "3-180/3-10" = c("3-10", "3-180")))

samp_order <- "station_depth"

propreads_bp <- ggplot(division_pivot_wmeta4, aes(x=reorder(get(.env$samp_order), desc(get(.env$samp_order))), y=value, fill = Taxonomic_group))+
  #labs(title = "Proportional read abundance")+
  theme_classic()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(month), cols = vars(newsizefract), scales = "free_y", space = "free_y")+
  ylab("Proportion of reads")+
  xlab("Station - Depth (m)")+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.1, 0, 0.1, "cm")))+
  #ylim(0,1.001)+
  ylim(0,1.01)+
  coord_flip()+
  theme(legend.title = element_text(size = 9),
                     legend.text = element_text(size = 8),
                     legend.key.size = unit(0.65, "lines"))+
  labs(fill = "Division/Class")+
  theme(legend.position = "bottom")
propreads_bp

annotate_figure(propreads_bp, top = text_grob(expression(paste("Size fraction (", mu, "m)"))), right = text_grob("Cruise month                 ", rot = -90))

