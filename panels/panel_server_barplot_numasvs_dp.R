#### Output ####
output$barplot_12 <- renderPlot({
  barplots()$p12
})

output$barplot_345 <- renderPlot({
  barplots()$p345
})

output$barplot_net <- renderPlot({
  barplots()$pnet
})


output$barplot_asvly_12 <- renderPlotly({
  barplots_asv()$p12ly
})

output$barplot_asvly_345 <- renderPlotly({
  barplots_asv()$p345ly
})

output$barplot_asvly_net <- renderPlotly({
  barplots_asv()$pnetly
})

barplots_asv <- eventReactive(input$actionb_barplot_asv, { 
  
  if (input$which_tab_rich == "merged") {
    asvtab_bp_start <- asvtab6_merged_subsamp_pa
  } else {
    asvtab_bp_start <- asvtab3_nonmerged_pa
  }
  
  #### Function to identify taxa that were present with < 5% of the reads in all samples ####
  limfun <- function(x) {
    ifelse(x>=input$prop_lim_asv/100,1,0)
  }
  
  taxlevel_asv <- input$taxlevel_barplot_asv
  
  #### Filter asv table to keep selected Divisions, unless 'All' is selected ####
  if (input$taxo_group_barp_asv != "All") {
    asvtab_filter <- asvtab_bp_start %>% filter(divisionlong %in% input$taxo_group_barp_asv)
  } else {
    asvtab_filter <- asvtab_bp_start
  }
  
  asvtab_filter_num <- asvtab_filter %>% purrr::keep(is.numeric)
  asvtab_filter_notnum <- asvtab_filter %>% purrr::keep(negate(is.numeric))
  
  #### If selected, create proportions ####
  if (input$propchoice_asv == "yes") {
    asvtab_filter_prop <- sweep(asvtab_filter_num, 2 , colSums(asvtab_filter_num), FUN = "/")
  } else {
    asvtab_filter_prop <- asvtab_filter_num
  }
  
  #### combine asv table with new proportions with taxonomy ####
  asvtab_filter_wtax <- cbind.data.frame(asvtab_filter_notnum, asvtab_filter_prop)
  
  asvtab_filter_wtax[is.na(asvtab_filter_wtax)] <- 0 #NA's because some samples have 0% of a certain class 
  
  #### Group by selected taxonomic level ####
  groupby_taxlevel <-  asvtab_filter_wtax %>% group_by(.data[[taxlevel_asv]]) %>% summarise_if(is.numeric, sum)
  
  ## Transform grouped df to matrix, to be able to apply limfun ##
  taxgroupspre_mat <- as.matrix(groupby_taxlevel[,-1, drop = FALSE])
  taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun)
  
  if (is.vector(taxgroupspre_bin)) {
    taxgroupspre_bin2 <- as.data.frame(as.list(taxgroupspre_bin))
  } else {
    taxgroupspre_bin2 <- as.data.frame(taxgroupspre_bin)}
  
  #### Identify taxa present as < selected % in all samples (low prop. taxa) ####
  taxgroupspre_bin_sums <- taxgroupspre_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
  taxgroupspre_bin_sums$taxgroups <- groupby_taxlevel %>% pull(get(taxlevel_asv))
  taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% select(taxgroups)
  
  taxgroups_select <- groupby_taxlevel %>% filter(get(taxlevel_asv) %in% taxgroupspre_bin_sums_yes$taxgroups)
  taxgroups_other <- groupby_taxlevel %>% filter(!get(taxlevel_asv) %in% taxgroupspre_bin_sums_yes$taxgroups)
  
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
    
    taxlevel_pivot_wmeta2_jan_mar <- taxlevel_pivot_wmeta2 %>% filter(month %in% c("Jan", "Mar"), fraction_max != 200)
    
    p12 <- ggplot(taxlevel_pivot_wmeta2_jan_mar, aes(x=reorder(station_depth, desc(station_depth)), y=value, fill = Taxonomic_group, 
                                                     text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_depth, percentornum, roundfun(value))))+
      labs(title = "Proportional read abundance")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    p12
    p12ly <- ggplotly(p12, tooltip = "text")
    
    #### Plot for May, Aug, Nov, Niskin samples ####
    taxlevel_pivot_wmeta2_may_nov <- taxlevel_pivot_wmeta2 %>% filter(month %in% c("May", "Aug", "Nov"), collection_method == "niskin")
    
    p345 <- ggplot(taxlevel_pivot_wmeta2_may_nov, aes(x=reorder(station_depth, desc(station_depth)), y=value, fill = Taxonomic_group,
                                                      text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_depth, percentornum, roundfun(value))))+
      #labs(title = "Proportional read abundance")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    
    p345ly <- ggplotly(p345, tooltip = "text")
    #### Plot net haul ####
    taxlevel_pivot_wmeta2_net <- taxlevel_pivot_wmeta2 %>% filter(collection_method == "net")
    
    pnet <- ggplot(taxlevel_pivot_wmeta2_net, aes(x=reorder(station_depth, desc(station_depth)), y=value, fill = Taxonomic_group,
                                                  text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_depth, percentornum, roundfun(value))))+
      #labs(title = "Proportional read abundance")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    pnetly <- ggplotly(pnet, tooltip = "text")
    
  } else {
    
    taxlevel_pivot_wmeta2 <- taxlevel_pivot_wmeta %>% mutate(station_dep_com = factor(station_dep_com, levels = unique(meta_tab_seq_event$station_dep_com), ordered = T),
                                                             month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov"), ordered = T), 
                                                             collection_method = factor(collection_method, levels = c("niskin", "net")),
                                                             size_fraction = factor(size_fraction, levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T))
    
    taxlevel_pivot_wmeta2_jan_mar <- taxlevel_pivot_wmeta2 %>% filter(month %in% c("Jan", "Mar"), fraction_max != 200)
    
    p12 <- ggplot(taxlevel_pivot_wmeta2_jan_mar, aes(x=reorder(station_dep_com, desc(station_dep_com)), y=value, fill = Taxonomic_group, 
                                                     text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_dep_com, percentornum, roundfun(value))))+
      labs(title = "Number or proportion of ASVs")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    p12
    p12ly <- ggplotly(p12, tooltip = "text")
    
    #### Plot for May, Aug, Nov, Niskin samples ####
    taxlevel_pivot_wmeta2_may_nov <- taxlevel_pivot_wmeta2 %>% filter(month %in% c("May", "Aug", "Nov"), collection_method == "niskin")
    
    p345 <- ggplot(taxlevel_pivot_wmeta2_may_nov, aes(x=reorder(station_dep_com, desc(station_dep_com)), y=value, fill = Taxonomic_group,
                                                      text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_dep_com, percentornum, roundfun(value))))+
      #labs(title = "Proportional read abundance")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    
    p345ly <- ggplotly(p345, tooltip = "text")
    #### Plot net haul ####
    taxlevel_pivot_wmeta2_net <- taxlevel_pivot_wmeta2 %>% filter(collection_method == "net")
    
    pnet <- ggplot(taxlevel_pivot_wmeta2_net, aes(x=reorder(station_dep_com, desc(station_dep_com)), y=value, fill = Taxonomic_group,
                                                  text = sprintf("Taxon: %s<br>Sample: %s<br>%s: %s ", Taxonomic_group, station_dep_com, percentornum, roundfun(value))))+
      #labs(title = "Proportional read abundance")+
      geom_bar(stat = "identity", position = "stack")+
      #scale_fill_manual(values = bpcolvec)+
      facet_grid(rows = vars(month), cols = vars(size_fraction), scales = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      #ylim(0,1.01)+
      coord_flip()
    pnetly <- ggplotly(pnet, tooltip = "text")
    
  }
  
  return(list(p12 = p12, p345 = p345, pnet = pnet, p12ly = p12ly, p345ly = p345ly, pnetly = pnetly))
})
