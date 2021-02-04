sbp_rare = sidebarPanel(
  radioButtons("which_tab", "ASV-table with separate or merged replicates:", 
               choices = c("Merged" = "merged", "Separate" = "separate"), selected = "merged"),
  selectInput("taxo_group_raref", label = "Division", 
                     choices = c(levels(as.factor(asvtab_merg_subsamp_readnum$divisionlong)),  "All"), 
                     selected = "All", multiple = TRUE),
  radioButtons("subsampled", "Subsampled to equal read number?", choices = c("Yes", "No"), selected = "Yes")
  
  
)

rarefpage = fluidPage(
  headerPanel("Rarefaction curves"),
  sidebarLayout(
    sbp_rare,
    mainPanel(
      plotOutput(outputId = "rarefactionplot_dp")
    )
  )
)
