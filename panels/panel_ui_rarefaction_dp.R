sbp_rare = sidebarPanel(
  actionButton("actionb_rarefcurve", "Create graphic", icon("signal")),
  radioButtons("which_asvtab", "ASV-table with separate or merged replicates:", 
               choices = c("Separate" = "separate", "Merged and subsampled" = "merged"), selected = "merged"),
  selectInput("taxo_group_raref", label = "Division", 
                     choices = c(levels(as.factor(asvtab3_merged_subsamp_readnum$divisionlong)),  "All"), 
                     selected = "All", multiple = TRUE)
  
  
)

rarefpage = fluidPage(
  waiter::use_waiter(),
  headerPanel("Rarefaction curves"),
  sidebarLayout(
    sbp_rare,
    mainPanel(
      plotlyOutput(outputId = "rarefactionplot_dp"),
      plotlyOutput(outputId = "nasvs_slope")
    )
  )
)
