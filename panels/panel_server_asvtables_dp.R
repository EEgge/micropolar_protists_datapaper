#### Output ####
output$tableS1 <- renderTable({
  asvtables()$tableS1
}, digits = 0)

asvtables <- reactive({
#### basic stats ####
taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")

# Total number of ASVs per Division #
nASVs_taxlevel <- asvtab5_merged_subsamp_prop %>% group_by(.data[[input$taxlevel_asvtables]]) %>% count()
# Number of ASVs per Division in each size fraction #
print(sum(nASVs_taxlevel[,2]))

asvtab_pa_longer <- pivot_longer(asvtab6_merged_subsamp_pa, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
asvtab_pa_longer_meta <- left_join(asvtab_pa_longer, meta_tab_sample_sf, by = "sample_sizefract")


numotus_asv_sizefract <- asvtab_pa_longer_meta %>% group_by(size_fraction, .data[[input$taxlevel_asvtables]], asv_code) %>% 
  summarise(sumval = sum(value)) %>% mutate(presabs = ifelse(sumval>0, 1, 0))
numotus_division_sizefract <- numotus_asv_sizefract %>% group_by(size_fraction, .data[[input$taxlevel_asvtables]]) %>% summarise(nasvs = sum(presabs))

roundfun <- function (x) {round(x,0)}

numotus_division_sizefract_wide <- pivot_wider(numotus_division_sizefract, id_cols = .data[[input$taxlevel_asvtables]], names_from = "size_fraction", values_from = "nasvs" )
numotus_division_sizefract_wtotal <- left_join(numotus_division_sizefract_wide, nASVs_taxlevel, by = .data[[input$taxlevel_asvtables]]) %>% arrange(desc(n)) %>% 
  select(.data[[input$taxlevel_asvtables]], '0.4_3', '3_10', '10_50', '50_200', '3_180', 'n') %>% 
  mutate_if(is.numeric, roundfun)
return(list(tableS1 = numotus_division_sizefract_wtotal))

})
#Min-max per division #
# asvtab_long <- pivot_longer(asvtab_subsamp_prop_tax, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
# asvtab_long_meta <- left_join(asvtab_long, meta_tab_nodups, by = "sample_sizefract")
# asvtab_sumdiv <- asvtab_long_meta %>% group_by(size_fraction, divisionlong, sample_sizefract) %>% summarise(prop_div = sum(value)) 
# asvtab_sumdiv_minmax <- asvtab_sumdiv %>% group_by(size_fraction, divisionlong) %>% summarise(minval = min(prop_div), maxval = max(prop_div)) %>% 
#   mutate(minmax = paste(round(minval,3), round(maxval,3), sep = ", "))
# asvtab_sumdiv_minmax_wide <- pivot_wider(asvtab_sumdiv_minmax, id_cols = "divisionlong", names_from = "size_fraction", values_from = "minmax") %>% 
#   left_join(nASVs_divisionlong, by ="divisionlong") %>% arrange(desc(n)) %>% select(divisionlong, '0.4_3', '3_10', '10_50', '50_200', '3_180', -n)
# 
# print(xtable(asvtab_sumdiv_minmax_wide, type = "latex", digits = 0), file = "propreads_division_sf.tex", include.rownames = FALSE)
