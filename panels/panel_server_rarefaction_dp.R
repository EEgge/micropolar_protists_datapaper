#### Output ####
output$rarefactionplot_dp <- renderPlotly({
  rarefaction_curves()$rareplotly
})

output$nasvs_slope <- renderPlotly({
  rarefaction_curves()$nasvs_slope_plotly
})


rarefaction_curves <- eventReactive(input$actionb_rarefcurve, {
  
  
  waiter <- waiter::Waiter$new()
  waiter$show()
  on.exit(waiter$hide())

  Sys.sleep(sample(5, 1))
  runif(1)

  if (input$which_asvtab == "separate") {
    asvtab_start <- asvtab1_nonmerged_readnum 
  } else {
    asvtab_start <- asvtab3_merged_subsamp_readnum
  }
  
  
  if (input$taxo_group_raref != "All") {
    asvtab_num <- asvtab_start %>% filter(divisionlong %in% input$taxo_group_raref) %>% 
      purrr::keep(is.numeric) #Create new ASV table with only the selected division
  } else {
    asvtab_num <- asvtab_start %>% purrr::keep(is.numeric)
  }
  
  rareobj <- vegan::rarecurve(t(asvtab_num),step=1000) 

rareo_msamp <- list()

if (input$which_asvtab == "separate") {
  for (i in 1:length(rareobj)) {
    rareo_msamp[[i]] <- cbind.data.frame("nASV" = rareobj[[i]], "step" = attr(rareobj[[i]], "Subsample"), "seq_event" =rep(names(asvtab_num)[i], length(rareobj[[i]])))
  }
  print(dim(rareo_msamp[[1]]))
  
  rareo_msamp_tab <- dplyr::bind_rows(rareo_msamp, .id = "column_label")
  
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

rareo_msamp_tab <- dplyr::bind_rows(rareo_msamp, .id = "column_label")

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
  rarslopes[i] <- vegan::rareslope(t(asvtab_num[,i]), sample=sum(asvtab_num[,i])-1)
}
names(rarslopes) <- names(asvtab_num)

min(rarslopes)

asvtab3c_num <- asvtab3c_merged_subsamp_pa %>% purrr::keep(is.numeric)
nasvs_sample_sizefract <- tibble("sample_sizefract" = names(asvtab3c_num), "nasvs" = colSums(asvtab3c_num))
slope_sample_sizefract <- tibble("sample_sizefract" = names(rarslopes), "slope" = rarslopes)

nasvs_slope <- left_join(nasvs_sample_sizefract, slope_sample_sizefract, by = "sample_sizefract") %>% rowwise() %>% 
  mutate(sizefract = paste(str_split(sample_sizefract, "_")[[1]][4],str_split(sample_sizefract, "_")[[1]][5], sep = "_"))

nasvs_slope_plot <- ggplot(data = nasvs_slope, aes(x = nasvs, y = slope, colour = sizefract, text = sprintf("Sample: %s<br>Number of ASVs: %s<br>Slope: %s ", sample_sizefract, round(nasvs,0), round(slope, 5))))+
  geom_point()

nasvs_slope_plotly <- ggplotly(nasvs_slope_plot, tooltip = "text")

  return(list(rareplotly = rareplotly, nasvs_slope_plotly = nasvs_slope_plotly))   
})
