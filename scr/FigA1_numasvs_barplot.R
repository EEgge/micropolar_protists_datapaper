asvtab6 <- read_delim(here("data", "asvtab6_merged_subsamp_pa.txt"), delim = "\t")
asvtab6_num <- asvtab6 %>% purrr::keep(is.numeric)
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]


asvtab_bp_start <- asvtab6_merged_subsamp_pa

taxgroups_abund <- c("Cercozoa", "Chlorophyta", "Choanoflagellida", "Ciliophora", "Cryptophyta", "Dinoflagellata_Dinophyceae", "Dinoflagellata_Syndiniales", "Haptophyta",
  "Katablepharidophyta", "Ochrophyta_Bacillariophyta", "Ochrophyta_Dictyochophyceae", "Ochrophyta_Pelagophyceae", "Picozoa", "Pseudofungi", "Radiolaria", "Sagenista")             

#### Group by selected taxonomic level ####
groupby_taxlevel <-  asvtab_filter_wtax %>% group_by(divisionlong) %>% summarise_if(is.numeric, sum)


#### Identify taxa present as < selected % in all samples (low prop. taxa) ####

taxgroups_select <- groupby_taxlevel %>% filter(divisionlong %in% taxgroups_abund)
taxgroups_other <- groupby_taxlevel %>% filter(!divisionlong %in% taxgroups_abund)

if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
  othersum <- colSums(taxgroups_other[,-1])
} else {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
  othersum <- NULL
}


taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)


#### Pivot table ####
if (input$which_tab_rich == "merged") {
  taxlevel_pivot <- pivot_longer(taxgroups_select3, cols = !Taxonomic_group, names_to = "sample_sizefract")
} else {
  taxlevel_pivot <- pivot_longer(taxgroups_select3, cols = !Taxonomic_group, names_to = "seq_event")
}

#### Join pivoted taxa table with metadata ####
if (input$which_tab_rich == "merged") {
  taxlevel_pivot_wmeta <- left_join(taxlevel_pivot, meta_tab_sample_sf, by = "sample_sizefract")
} else {
  taxlevel_pivot_wmeta <- left_join(taxlevel_pivot, meta_tab_seq_event, by = "seq_event")
}

### 8.2 her er jeg
if (input$propchoice_asv == "yes") {
  roundfun <- function(x) {round(x*100, 3)}
} else {
  roundfun <- function(x) {round(x,0)}
}

if (input$propchoice_asv == "yes") {
  percentornum <- "Percent"
} else {
  percentornum <- "Number of ASVs"
}


#### 03.02.21 her er jeg ####
if (input$which_tab_rich == "merged") {
  taxlevel_pivot_wmeta2 <- taxlevel_pivot_wmeta %>% mutate(station_depth = factor(station_depth, levels = unique(meta_tab_sample_sf$station_depth), ordered = T),
                                                           month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T), 
                                                           collection_method = factor(collection_method, levels = c("niskin", "net")), 
                                                           size_fraction = factor(size_fraction, levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T))
  taxlevel_pivot_wmeta3 <- taxlevel_pivot_wmeta2 %>% mutate(size_fraction=recode(size_fraction, "0.4_3" = "0.45-3", "3_10" = "3-10", "10_50" = "10-50", "50_200" = "50-200", "3_180" = "3-180"))
  taxlevel_pivot_wmeta4 <- taxlevel_pivot_wmeta3 %>% mutate(newsizefract = forcats::fct_collapse(size_fraction, "3-180/3-10" = c("3-10", "3-180")))
  
  taxlevel_pivot_wmeta2_jan_mar <- taxlevel_pivot_wmeta2 %>% filter(month %in% c("Jan", "Mar"), fraction_max != 200)
  
  p12 <- ggplot(taxlevel_pivot_wmeta4, aes(x=reorder(station_depth, desc(station_depth)), y=value, fill = Taxonomic_group, 
                                           text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_depth, percentornum, roundfun(value))))+
    labs(title = percentornum)+
    geom_bar(stat = "identity", position = "stack")+
    #scale_fill_manual(values = bpcolvec)+
    facet_grid(rows = vars(month), cols = vars(newsizefract), scales = "free_y", space = "free_y")+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    #ylim(0,1.01)+
    coord_flip()
  p12
  