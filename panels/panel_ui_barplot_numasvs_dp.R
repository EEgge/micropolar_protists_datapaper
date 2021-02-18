sbp_barasvs = sidebarPanel(
  actionButton("actionb_barplot_asv", "Create graphic", icon("signal")),
  numericInput("prop_lim_asv", "Taxonomic groups present with < x% in all samples will be grouped together as 'Other'",
               value = 5, min = 0, max = 100),
  radioButtons("which_tab_rich", "ASV-table with merged or separate replicates:", 
               choices = c("Merged" = "merged", "Separate" = "separate"), selected = "merged"),
  selectInput("taxo_group_barp_asv", label = "Division", 
              choices = c(levels(as.factor(asvtab4_merged_subsamp_readnum$divisionlong)),  "All"), 
              selected = "All", multiple = TRUE),
  radioButtons("taxlevel_barplot_asv", label = "Taxonomic level",
               choices = c("Division" = "divisionlong", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species")),
  radioButtons("propchoice_asv", label = "Proportion", choices = c("Yes" = "yes","No" = "no"))
  
)


barpage_asvs = fluidPage(
  waiter::use_waiter(),
  headerPanel("Barplot, taxonomic distribution of ASV richness"),
  sidebarLayout(
    sbp_barasvs,
    mainPanel(
      plotlyOutput(outputId = "barplot_asvly_12"),
      plotlyOutput(outputId = "barplot_asvly_345"),
      plotlyOutput(outputId = "barplot_asvly_net")
    )
  )
)

