#### Output ####
output$rarefactionplot_dp <- renderPlotly({
  rarefaction_curves()$rareplotly
})

rarefaction_curves <- eventReactive(input$actionb_rarefcurve, {
  # waiter::Waiter$new(id = "rarefactionplot_dp")$show()
  # 
  # Sys.sleep(3)
  # 
  
  
  
  waiter <- waiter::Waiter$new()
  waiter$show()
  on.exit(waiter$hide())

  Sys.sleep(sample(5, 1))
  runif(1)

  if (input$which_asvtab == "separate") {
    asvtab_start <- asvtab1_nonmerged_readnum 
  } else {
    asvtab_start <- asvtab4_merged_subsamp_readnum
  }
  
  
  if (input$taxo_group_raref != "All") {
    asvtab_num <- asvtab_start %>% filter(divisionlong %in% input$taxo_group_raref) %>% 
      purrr::keep(is.numeric) #Create new OTU table with only the selected division
  } else {
    asvtab_num <- asvtab_start %>% purrr::keep(is.numeric)
  }
  
  rareobj <- rarecurve(t(asvtab_num),step=1000) 

rareo_msamp <- list()

if (input$which_asvtab == "separate") {
  for (i in 1:length(rareobj)) {
    rareo_msamp[[i]] <- cbind.data.frame("nASV" = rareobj[[i]], "step" = attr(rareobj[[i]], "Subsample"), "seq_event" =rep(names(asvtab_num)[i], length(rareobj[[i]])))
  }
  print(dim(rareo_msamp[[1]]))
  
  rareo_msamp_tab <- bind_rows(rareo_msamp, .id = "column_label")
  
  rareo_msamp_tab_meta <- left_join(rareo_msamp_tab, meta_tab_seq_event, by = "seq_event")
  
  rareo_msamp_tab_meta <- rareo_msamp_tab_meta %>% mutate_at(vars(size_fraction), list( ~factor(., levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T)))
  
  rareplot <- ggplot(rareo_msamp_tab_meta, aes(x = step, y = nASV, group = seq_event, 
                                               text = sprintf("Sample: %s<br>Number of ASVs: %s ", seq_event, round(nASV,0))))+
    geom_line(aes(color=size_fraction))+
    xlab("Number of reads")+
    ylab("Number of ASVs")
  
  rareplotly <- ggplotly(rareplot, tooltip = "text")
  
  
} else {
for (i in 1:length(rareobj)) {
rareo_msamp[[i]] <- cbind.data.frame("nASV" = rareobj[[i]], "step" = attr(rareobj[[i]], "Subsample"), "sample_sizefract" =rep(names(asvtab_num)[i], length(rareobj[[i]])))
}
print(dim(rareo_msamp[[1]]))

rareo_msamp_tab <- bind_rows(rareo_msamp, .id = "column_label")

rareo_msamp_tab_meta <- left_join(rareo_msamp_tab, meta_tab_sample_sf, by = "sample_sizefract")

rareo_msamp_tab_meta <- rareo_msamp_tab_meta %>% mutate_at(vars(size_fraction), list( ~factor(., levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T)))

rareplot <- ggplot(rareo_msamp_tab_meta, aes(x = step, y = nASV, group = sample_sizefract, 
                                             text = sprintf("Sample: %s<br>Number of ASVs: %s ", sample_sizefract, round(nASV,0))))+
  geom_line(aes(color=size_fraction))+
  xlab("Number of reads")+
  ylab("Number of ASVs")

rareplotly <- ggplotly(rareplot, tooltip = "text")

}


rarslopes <- c(NULL)
for (i in 1:dim(asvtab_num)[2]) {
  rarslopes[i] <- rareslope(t(asvtab_num[,i]), sample=sum(asvtab_num[,i])-1)
}
names(rarslopes) <- names(asvtab_num)

min(rarslopes)

asvtab6_num <- asvtab6_merged_subsamp_pa %>% purrr::keep(is.numeric)
nasvs_sample_sizefract <- tibble("sample_sizefract" = names(asvtab6_num), "nasvs" = colSums(asvtab6_num))
slope_sample_sizefract <- tibble("sample_sizefract" = names(rarslopes), "slope" = rarslopes)

nasvs_slope <- left_join(nasvs_sample_sizefract, slope_sample_sizefract, by = "sample_sizefract")

  return(list(rareplotly = rareplotly))   
})
