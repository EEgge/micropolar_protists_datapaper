require(xtable)
taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "class", "family", "order", "genus", "species", "sequence", "divisionlong")
asvtab <- read_delim(here("data", "asvtab6_merged_subsamp_pa.txt"), delim = "\t")
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles_forplot.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

nASVs_divisionlong <- asvtab %>% group_by(divisionlong) %>% count()

# pivot and add meta data to presence-absence #
asvtab_pa_longer <- pivot_longer(asvtab, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
asvtab_pa_longer_meta <- left_join(asvtab_pa_longer, meta_tab_sample_sf, by = "sample_sizefract")

# Number of ASVs per Division in each size fraction #
numotus_asv_sizefract <- asvtab_pa_longer_meta %>% group_by(size_fraction, divisionlong, asv_code) %>% 
  summarise(sumval = sum(value)) %>% mutate(presabs = ifelse(sumval>0, 1, 0))
numotus_division_sizefract <- numotus_asv_sizefract %>% group_by(size_fraction, divisionlong) %>% summarise(nasvs = sum(presabs))

numotus_division_sizefract_wide <- pivot_wider(numotus_division_sizefract, id_cols = "divisionlong", names_from = "size_fraction", values_from = "nasvs" )
numotus_division_sizefract_wtotal <- left_join(numotus_division_sizefract_wide, nASVs_divisionlong, by = "divisionlong") %>% arrange(desc(n)) %>% 
  select(divisionlong, '0.4_3', '3_10', '10_50', '50_200', '3_180', 'n')

print(xtable(numotus_division_sizefract_wtotal, type = "latex", digits = 0), file = "figures_tables/tables/TableA1_numasvs_division_sf.tex", include.rownames = FALSE)
