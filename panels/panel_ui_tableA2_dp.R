sbp_tableS2 = sidebarPanel(
  radioButtons("taxlevel_tableS2", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"))
  
)

tableS2page = fluidPage(
  headerPanel("Table A2 - Min. and max. percentage of reads assigned to each taxonomic group, distributed by size fraction. The entries have the format 'min., max.'"),
  sidebarLayout(
    sbp_tableS2,
    mainPanel(
      DT::dataTableOutput(outputId = "tableS2")
    )
  )
)
