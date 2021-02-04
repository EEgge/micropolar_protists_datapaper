#### Output ####
output$rarefactionplot_dp <- renderPlot({
  rarefaction_curves()$rareplot
})

rarefaction_curves <- reactive ({
  if (input$taxo_group_raref != "All") {
    asvtab_num <- asvtab_merg_subsamp_readnum %>% filter(divisionlong %in% input$taxo_group_raref) %>% 
      purrr::keep(is.numeric) #Create new OTU table with only the selected division
  } else {
    asvtab_num <- asvtab_merg_subsamp_readnum %>% purrr::keep(is.numeric)
  }
  
  rareobj <- rarecurve(t(asvtab_num),step=1000) 

rareo_msamp <- list()
for (i in 1:length(rareobj)) {
rareo_msamp[[i]] <- cbind.data.frame("nASV" = rareobj[[i]], "step" = attr(rareobj[[i]], "Subsample"), "sample_sizefract" =rep(names(asvtab_num)[i], length(rareobj[[i]])))
}
print(dim(rareo_msamp[[1]]))

rareo_msamp_tab <- bind_rows(rareo_msamp, .id = "column_label")

rareo_msamp_tab_meta <- left_join(rareo_msamp_tab, meta_tab_sample_sf, by = "sample_sizefract")

rareo_msamp_tab_meta <- rareo_msamp_tab_meta %>% mutate_at(vars(size_fraction), list( ~factor(., levels = c("0.4_3", "3_10", "10_50", "50_200", "3_180"), ordered = T)))

rareplot <- ggplot(rareo_msamp_tab_meta, aes(x = step, y = nASV, group = sample_sizefract))+
  geom_line(aes(color=size_fraction))+
  xlab("Number of reads")+
  ylab("Number of ASVs")

rarslopes <- c(NULL)
for (i in 1:dim(asvtab_num)[2]) {
  rarslopes[i] <- rareslope(t(asvtab_num[,i]), sample=sum(asvtab_num[,i])-1)  
}
names(rarslopes) <- names(asvtab_num)

min(rarslopes)

nasvs_sample_sizefract <- tibble("sample_sizefract" = names(asvtab_subsamp_num_prop_pa), "nasvs" = colSums(asvtab_subsamp_num_prop_pa))
slope_sample_sizefract <- tibble("sample_sizefract" = names(rarslopes), "slope" = rarslopes)

nasvs_slope <- left_join(nasvs_sample_sizefract, slope_sample_sizefract, by = "sample_sizefract")

  return(list(rareplot = rareplot))   
})
