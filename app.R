#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(vegan)
library(here)
library(plotly)
library(ggpubr)
library(data.table)
library(vegan)
library(readr)


#### Read files ####
## ASV tables ##
asvtab1_nonmerged_readnum <- read_delim(here("data", "asvtab1_nonmerged_readnum.txt"), delim = "\t")
asvtab2_nonmerged_prop <- read_delim(here("data", "asvtab2_nonmerged_prop.txt"), delim = "\t")
asvtab3_nonmerged_pa <- read_delim(here("data", "asvtab3_nonmerged_pa.txt"), delim = "\t")
asvtab4_merged_subsamp_readnum <- read_delim(here("data", "asvtab4_merged_subsamp_readnum.txt"), delim = "\t")
asvtab5_merged_subsamp_prop <- read_delim(here("data", "asvtab5_merged_subsamp_prop.txt"), delim = "\t")
asvtab6_merged_subsamp_pa <- read_delim(here("data", "asvtab6_merged_subsamp_pa.txt"), delim = "\t")

## meta data for each sequencing_event ##
meta_tab_seq_event <- read_delim(here("data", "meta_data_fastqfiles.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab_seq_event$sample_sizefract),]

## Environmental data files ##
env_table <- read_delim(here("data", "env_data_depths.txt"), delim = "\t")
env_profiles <- read_delim(here("data", "env_data_profiles.txt"), delim = "\t")


source("panels/panel_ui_barplot_dp.R", local = TRUE)
source("panels/panel_ui_rarefaction_dp.R", local = TRUE)
source("panels/panel_ui_barplot_numasvs_dp.R", local = TRUE)
source("panels/panel_ui_env_data_dp.R", local = TRUE)
source("panels/panel_ui_asvtables_dp.R", local = TRUE)
source("panels/panel_ui_tableS2_dp.R", local = TRUE)

ui <- navbarPage(
    "Interactive figures from MicroPolar protist metabarcoding data paper",
    tabPanel("Overview", includeMarkdown("panels/description_figures_tables.Rmd")),
    #tabPanel("Figure 1 - Map", mappage)
    tabPanel("Figure 2 - Environmental data", envpage),
    tabPanel("Figure 3 - Rarefaction curves", rarefpage),
    tabPanel("Figure 4 - Barplot - proportion of reads", barpage),
    tabPanel("Figure S1 - Barplot - number of ASVs", barpage_asvs),
    tabPanel("Table S1 - Number of ASVs per clade in each size fraction", asvtablespage),
    tabPanel("Table S2 - Min. and max. proportion of reads per clade in each size fraction", tableS2page)
    
)

# Define server logic
server <- function(input, output) {
    
    #Environmental data
    source("panels/panel_server_env_data_dp.R", local = TRUE)
    #Rarefaction panel
    source("panels/panel_server_rarefaction_dp.R", local = TRUE)
    #Barplot panel
    source("panels/panel_server_barplot_dp.R", local = TRUE)
    source("panels/panel_server_barplot_numasvs_dp.R", local = TRUE)
    
    source("panels/panel_server_asvtables_dp.R", local = TRUE)
    source("panels/panel_server_tableS2_dp.R", local = T)
    
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
