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
asvtab_wtax <- read_excel(here("metapr2/export", "metapr2_wide_asv_set_207_208_209_Eukaryota.xlsx"))
asvtab_subsamp_prop_tax <- read_delim(here("data", "asvtab_subsamp_prop_wtax.txt"), delim = "\t") #Metazoa and Embryophyceae removed, seq_event replicate sample_sf merged, subsampled to equal read number within each size fraction.

taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")

n_occur <- data.frame(table(asvtab_wtax$asv_code))
n_occur %>% filter(Freq > 1)
asvtab_wtax %>% filter(asv_code == "41a8ef8f8066c541e9cd9d3df7734166370a68fa") %>% purrr::keep(is.numeric) %>% sum()


asvtab_wtax <- asvtab_wtax %>% add_column(divisionlong = NA)
for (i in 1:dim(asvtab_wtax)[1]) {
  if (asvtab_wtax$division[i] %in% c("Dinoflagellata", "Ochrophyta")) {
    asvtab_wtax$divisionlong[i] <- paste(asvtab_wtax$division[i], asvtab_wtax$class[i], sep = "_")
  } else {
    asvtab_wtax$divisionlong[i] <- asvtab_wtax$division[i]
  }
}


taxlevels <- c(names(asvtab_wtax)[1:10], "divisionlong")
asvtab_wo_metazoa_embr <- asvtab_wtax %>% filter(division !="Metazoa", class != "Embryophyceae")
dim(asvtab_wo_metazoa_embr) #6536 x 209

#### Proportions of non-subsampled asv table ####
asvtab_num <- asvtab_wo_metazoa_embr %>% select_if(is.numeric)
asvtab_notnum <- asvtab_wo_metazoa_embr %>% select_if(negate(is.numeric))
asvtab_num_prop <- sweep(asvtab_num, 2 , colSums(asvtab_num), FUN = "/")

colSums(asvtab_num)
asvtab_prop_tax <- cbind.data.frame(asvtab_notnum, asvtab_num_prop)
#write_delim(asvtab_prop_tax, "data/asvtab_nonmerged_prop_wtax.txt", delim = "\t")





asvtab_subsamp_prop_tax <- read_delim(here("data", "asvtab5_merged_subsamp_prop.txt"), delim = "\t")
dim(asvtab_subsamp_prop_tax) #6536 x 166 

#### Create presence-absence after subsampling asv table ####
asvtab_subsamp_notnum <- asvtab_subsamp_prop_tax %>% select_if(negate(is.numeric))
asvtab_subsamp_num_prop <- asvtab_subsamp_prop_tax %>% select_if(is.numeric)
asvtab_subsamp_num_prop_pa <- asvtab_subsamp_num_prop
asvtab_subsamp_num_prop_pa[asvtab_subsamp_num_prop_pa >0] <- 1
asvtab_subsamp_pa_tax <- bind_cols(asvtab_subsamp_notnum, asvtab_subsamp_num_prop_pa)
#write_delim(asvtab_subsamp_pa_tax, "data\asvtab6_merged_subsamp_pa.txt", delim = "\t")

#Num ASVs after subsampling : 
length(which(rowSums(asvtab_subsamp_num_prop_pa) == 0))

# pivot and add meta data to presence-absence #
asvtab_pa_longer <- pivot_longer(asvtab_subsamp_pa_tax, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
asvtab_pa_longer_meta <- left_join(asvtab_pa_longer, meta_tab_sample_sf, by = "sample_sizefract")

#### basic stats ####
# Total number of ASVs per Division #
nASVs_divisionlong <- asvtab_subsamp_prop_tax %>% group_by(divisionlong) %>% count()
# Number of ASVs per Division in each size fraction #
numotus_asv_sizefract <- asvtab_pa_longer_meta %>% group_by(size_fraction, divisionlong, asv_code) %>% 
  summarise(sumval = sum(value)) %>% mutate(presabs = ifelse(sumval>0, 1, 0))
numotus_division_sizefract <- numotus_asv_sizefract %>% group_by(size_fraction, divisionlong) %>% summarise(nasvs = sum(presabs))

numotus_division_sizefract_wide <- pivot_wider(numotus_division_sizefract, id_cols = "divisionlong", names_from = "size_fraction", values_from = "nasvs" )
numotus_division_sizefract_wtotal <- left_join(numotus_division_sizefract_wide, nASVs_divisionlong, by = "divisionlong") %>% arrange(desc(n)) %>% 
  select(divisionlong, '0.4_3', '3_10', '10_50', '50_200', '3_180', 'n')

print(xtable(numotus_division_sizefract_wtotal, type = "latex", digits = 0), file = "numasvs_division_sf.tex", include.rownames = FALSE)

#Min-max per division #
asvtab_long <- pivot_longer(asvtab_subsamp_prop_tax, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
asvtab_long_meta <- left_join(asvtab_long, meta_tab_nodups, by = "sample_sizefract")
asvtab_sumdiv <- asvtab_long_meta %>% group_by(size_fraction, divisionlong, sample_sizefract) %>% summarise(prop_div = sum(value)) 
asvtab_sumdiv_minmax <- asvtab_sumdiv %>% group_by(size_fraction, divisionlong) %>% summarise(minval = min(prop_div), maxval = max(prop_div)) %>% 
  mutate(minmax = paste(round(minval,3), round(maxval,3), sep = ", "))
asvtab_sumdiv_minmax_wide <- pivot_wider(asvtab_sumdiv_minmax, id_cols = "divisionlong", names_from = "size_fraction", values_from = "minmax") %>% 
  left_join(nASVs_divisionlong, by ="divisionlong") %>% arrange(desc(n)) %>% select(divisionlong, '0.4_3', '3_10', '10_50', '50_200', '3_180', -n)
  
print(xtable(asvtab_sumdiv_minmax_wide, type = "latex", digits = 0), file = "propreads_division_sf.tex", include.rownames = FALSE)

##### BARPLOT #####
#### Make legend ####
bpcol <- read_delim(here("data","col_figS2.txt"), col_names = T, delim = "\t", comment = "")
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
groupby_division <- asvtab_subsamp_prop_tax %>% group_by(divisionlong) %>% summarise_if(is.numeric, sum)




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
  groupby_bin2 <- as.data.frame(groupby_bin)}


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

#### Create data frame with >= 5% taxa, and all < 5% taxa grouped together ##-##
taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)

#### Pivot table ####
division_pivot <- pivot_longer(taxgroups_select3, cols = !Taxonomic_group, names_to = "sample_sizefract")
#division_pivot <- pivot_longer(taxgroups_select3, cols = !Taxonomic_group, names_to = "seq_event")

## Transform Taxonomic_group to factor, change order of levels so "Eukaryota_unclassified" comes at the end ##
#Taxonomic_group2 <- Taxonomic_group[c(1:7,9:18,8,19)]
division_pivot <- division_pivot %>% mutate(Taxonomic_group = factor(Taxonomic_group, levels=.env$Taxonomic_group, ordered = T))

meta_tab_nodups <- meta_tab[!duplicated(meta_tab$sample_sizefract),]

#### Join pivoted taxa table with metadata ####
division_pivot_wmeta <- left_join(division_pivot, meta_tab_sample_sf, by = "sample_sizefract")
#division_pivot_wmeta <- left_join(division_pivot, meta_tab, by = "seq_event")

division_pivot_wmeta_merge <- division_pivot_wmeta %>% group_by(Taxonomic_group, sample_sizefract, station_depth, env_sample, size_fraction, month, station, collection_method) %>%
                              summarise_if(is.numeric, mean) %>% 
                              select(-dna_concentration, -`260_280`, -`260_230`, -sample_sequence)

division_pivot_wmeta2 <- division_pivot_wmeta %>% mutate(station_dep_com = factor(station_dep_com, levels = unique(meta_tab$station_dep_com), ordered = T),
                                                         month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T), 
                                                         collection_method = factor(collection_method, levels = c("niskin", "net")),
                                                         size_fraction = factor(size_fraction, levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T))

division_pivot_wmeta2 <- division_pivot_wmeta_merge %>% mutate(station_depth = factor(station_depth, levels = unique(meta_tab_sample_sf$station_depth), ordered = T),
                                                          month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T), 
                                                          collection_method = factor(collection_method, levels = c("niskin", "net")),
                                                          size_fraction = factor(size_fraction, levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T))

division_pivot_wmeta2 %>% ungroup %>% group_by(Taxonomic_group, size_fraction) %>% summarise_if(is.numeric)
#### Plot for January and March (only two size fractions, no net haul) ####
division_pivot_wmeta2_jan_mar <- division_pivot_wmeta2 %>% filter(month %in% c("Jan", "Mar"), fraction_max != 200)

samp_order <- "station_depth"

p12 <- ggplot(division_pivot_wmeta2_jan_mar, aes(x=reorder(get(.env$samp_order), desc(get(.env$samp_order))), y=value, fill = Taxonomic_group))+
  labs(title = "Proportional read abundance")+
  theme_classic()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
  xlab("")+
  ylab("")+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.1, 0, 0.1, "cm")))+
  ylim(0,1.001)+
  coord_flip()
p12

#### Plot for May, Aug, Nov, Niskin samples ####
division_pivot_wmeta2_may_nov <- division_pivot_wmeta2 %>% filter(month %in% c("May", "Aug", "Nov"), collection_method == "niskin")

p345 <- ggplot(division_pivot_wmeta2_may_nov, aes(x=reorder(get(.env$samp_order), desc(get(.env$samp_order))), y=value, fill = Taxonomic_group))+
  #labs(title = "Proportional read abundance")+
  theme_classic()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
  xlab("Station - Depth (m)")+
  ylab("")+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.1, 0, 0.1, "cm")))+
  ylim(0,1.001)+
  coord_flip()
p345

#### Plot net haul ####
division_pivot_wmeta2_net <- division_pivot_wmeta2 %>% filter(collection_method == "net")

pnet <- ggplot(division_pivot_wmeta2_net, aes(x=reorder(get(.env$samp_order), desc(get(.env$samp_order))), y=value, fill = Taxonomic_group))+
  #labs(title = "Proportional read abundance")+
  theme_classic()+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
  xlab("")+
  ylab("Proportion of reads")+
  theme(strip.background = element_rect(size = 0.5),
        strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")),
        strip.text.y = element_text(margin = margin(0, 0.1, 0, 0.1, "cm")))+
  ylim(0,1.001)+
  coord_flip()
pnet




# p12 <- ggplot(taxgroups_select3tdf_mp12melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+

#### Combine figures ####
figureS2 <- ggdraw()+
  draw_plot(p12+theme(legend.position = "none"), x = 0, y = 33/52, width = 0.545, height = 19/52)+
  draw_plot(p345+theme(legend.position = "none"), x = 0, y = 9.8/52, width = 1, height = 25/52)+
  draw_plot(pnet+theme(legend.position = "none"), x = 0.48, y = 0, width = 0.535, height = 11.5/52)#+
  #draw_grob(bp_legnd, x = 0.275, y = -.37)
figureS2
