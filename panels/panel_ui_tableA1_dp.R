sbp_asvtables = sidebarPanel(
  radioButtons("taxlevel_asvtables", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"))
  
)

asvtablespage = fluidPage(
  headerPanel("Table A1 - Number of ASVs assigned to each taxonomic group, distributed by size fraction."),
  sidebarLayout(
    sbp_asvtables,
    mainPanel(
      DT::dataTableOutput(outputId = "tableS1")
    )
  )
)
