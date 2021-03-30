sbp_bar = sidebarPanel(
    actionButton("actionb_barplot", "Create graphic", icon("signal")),
    numericInput("prop_lim", "Taxonomic groups present with < x% in all samples will be grouped together as 'Other'",
                 value = 5, min = 0, max = 100),
    radioButtons("which_tab", "ASV-table with merged or separate replicates:", 
                 choices = c("Merged" = "merged", "Separate" = "separate"), selected = "merged"),
    selectInput("taxo_group_barp", label = "Division", 
                choices = c(levels(as.factor(asvtab3_merged_subsamp_readnum$divisionlong)),  "All"), 
                selected = "All", multiple = TRUE),
    radioButtons("taxlevel_barplot", label = "Taxonomic level",
                 choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species")),
    radioButtons("propchoice", label = "Proportion of selected group or all", choices = c("Selected group" = "selectgroup","All" = "All"))
    
  )


barpage = fluidPage(
  waiter::use_waiter(),
  headerPanel("Barplot"),
  sidebarLayout(
    sbp_bar,
    mainPanel(
      plotlyOutput(outputId = "barplotly", height = "1200px", width = "1000px"),
      textOutput(outputId = "bp_text")
      
         )
  )
  )

