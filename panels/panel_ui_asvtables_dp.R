sbp_asvtables = sidebarPanel(
  # selectInput("taxo_group_barp_asv", label = "Division", 
  #             choices = c(levels(as.factor(asvtab4_merged_subsamp_readnum$divisionlong)),  "All"), 
  #             selected = "All", multiple = TRUE),
  radioButtons("taxlevel_asvtables", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"))
  
)

asvtablespage = fluidPage(
  headerPanel("Table "),
  sidebarLayout(
    sbp_asvtables,
    mainPanel(
      tableOutput(outputId = "tableS1")
    )
  )
)
