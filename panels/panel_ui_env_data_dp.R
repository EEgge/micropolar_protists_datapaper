envpage = fluidPage(
  headerPanel("Profiles of physical and chemical parameters"),
    mainPanel(
      plotlyOutput(outputId = "salplotly"),
      plotlyOutput(outputId = "tempplotly"),
      plotlyOutput(outputId = "sigmaplotly"),
      plotlyOutput(outputId = "inorgNplotly"),
      plotlyOutput(outputId = "phosplotly"),
      plotlyOutput(outputId = "silplotly"),
      plotlyOutput(outputId = "chlplotly")
      
    )
  )
