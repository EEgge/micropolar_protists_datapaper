sbp_tableS2 = sidebarPanel(
  # selectInput("taxo_group_barp_asv", label = "Division", 
  #             choices = c(levels(as.factor(asvtab4_merged_subsamp_readnum$divisionlong)),  "All"), 
  #             selected = "All", multiple = TRUE),
  radioButtons("taxlevel_tableS2", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"))
  
)

tableS2page = fluidPage(
  headerPanel("Table S2"),
  sidebarLayout(
    sbp_tableS2,
    mainPanel(
      tableOutput(outputId = "tableS2")
    )
  )
)
