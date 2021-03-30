sbp_tableS2 = sidebarPanel(
  radioButtons("taxlevel_tableS2", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"))
  
)

tableS2page = fluidPage(
  headerPanel("Table S2"),
  sidebarLayout(
    sbp_tableS2,
    mainPanel(
      DT::dataTableOutput(outputId = "tableS2")
    )
  )
)
