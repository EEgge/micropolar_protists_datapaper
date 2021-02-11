#### Output ####
output$tableS2 <- renderTable({
  minmaxtable()$tableS2
}, digits = 0)

minmaxtable <- reactive({
  #### basic stats ####
  taxlevels <- c("asv_code", "kingdom", "supergroup", "division", "divisionlong", "class", "family", "order", "genus", "species", "sequence")
  
  # Total number of ASVs per Division #
  nASVs_taxlevel <- asvtab5_merged_subsamp_prop %>% group_by(.data[[input$taxlevel_tableS2]]) %>% count()
  
  roundfun <- function (x) {round(x,0)}
  
  #Min-max per division #
  asvtabprop_long <- pivot_longer(asvtab5_merged_subsamp_prop, cols =!all_of(taxlevels),  names_to = "sample_sizefract")
  asvtabprop_long_meta <- left_join(asvtabprop_long, meta_tab_sample_sf, by = "sample_sizefract")
  asvtab_sumdiv <- asvtabprop_long_meta %>% group_by(size_fraction, .data[[input$taxlevel_tableS2]], sample_sizefract) %>% summarise(prop_div = sum(value))
  asvtab_sumdiv_minmax <- asvtab_sumdiv %>% group_by(size_fraction, .data[[input$taxlevel_tableS2]]) %>% summarise(minval = min(prop_div), maxval = max(prop_div)) %>%
    mutate(minmax = paste(round(minval,3), round(maxval,3), sep = ", "))
  asvtab_sum_minmax_wide <- pivot_wider(asvtab_sumdiv_minmax, id_cols = .data[[input$taxlevel_tableS2]], names_from = "size_fraction", values_from = "minmax") %>%
    left_join(nASVs_taxlevel, by =.data[[input$taxlevel_tableS2]]) %>% arrange(desc(n)) %>% select(.data[[input$taxlevel_tableS2]], '0.4_3', '3_10', '10_50', '50_200', '3_180', -n)
  
 
  return(list(tableS2 = asvtab_sum_minmax_wide))
  
})

