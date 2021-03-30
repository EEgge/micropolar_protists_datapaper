#Min-max per division #
require(xtable)
taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "class", "family", "order", "genus", "species", "sequence", "divisionlong")
asvtab <- read_delim(here("data", "asvtab5_merged_subsamp_prop.txt"), delim = "\t")
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

nASVs_divisionlong <- asvtab %>% group_by(divisionlong) %>% count()

asvtab_long <- pivot_longer(asvtab, cols =!all_of(taxlevels),  names_to = "sample_sizefract")

asvtab_long_meta <- left_join(asvtab_long, meta_tab_sample_sf, by = "sample_sizefract")
asvtab_sumdiv <- asvtab_long_meta %>% group_by(size_fraction, divisionlong, sample_sizefract) %>% summarise(prop_div = sum(value)) 
asvtab_sumdiv_minmax <- asvtab_sumdiv %>% group_by(size_fraction, divisionlong) %>% summarise(minval = min(prop_div), maxval = max(prop_div)) %>% 
mutate(minmax = paste(100*round(minval,3), 100*round(maxval,3), sep = ", "))
asvtab_sumdiv_minmax_wide <- pivot_wider(asvtab_sumdiv_minmax, id_cols = "divisionlong", names_from = "size_fraction", values_from = "minmax") %>% 
left_join(nASVs_divisionlong, by ="divisionlong") %>% arrange(desc(n)) %>% select(divisionlong, '0.4_3', '3_10', '10_50', '50_200', '3_180', -n)

print(xtable(asvtab_sumdiv_minmax_wide, type = "latex", digits = 0), file = "figures_tables/tables/TableA2_propreads_division_sf.tex", include.rownames = FALSE)
