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



# asvtab_wtax <- read_excel(here("metapr2/export", "metapr2_wide_asv_set_207_208_209_Eukaryota_nodups.xlsx"))
# asvtab_wtax <- asvtab_wtax %>% add_column(divisionlong = NA)
# for (i in 1:dim(asvtab_wtax)[1]) {
#     if (asvtab_wtax$division[i] %in% c("Dinoflagellata", "Ochrophyta")) {
#         asvtab_wtax$divisionlong[i] <- paste(asvtab_wtax$division[i], asvtab_wtax$class[i], sep = "_")
#     } else {
#         asvtab_wtax$divisionlong[i] <- asvtab_wtax$division[i]
#     }
# }

asvtab_merg_subsamp_readnum <- read_delim(here("data", "asv_tab_merged_subsamp.txt"), delim = "\t")
asvtab_subsamp_prop_tax <- read_delim(here("data", "asvtab_subsamp_prop_tax.txt"), delim = "\t")
meta_tab_seq_event <- read_delim(here("figures_tables/tables", "table1.txt"), delim = "\t")
meta_tab_sample_sf <- meta_tab_seq_event[!duplicated(meta_tab$sample_sizefract),]


source("panels/panel_ui_barplot_dp.R", local = TRUE)
source("panels/panel_ui_rarefaction_dp.R", local = TRUE)

ui <- navbarPage(
    "Interactive figures from MicroPolar protist metabarcoding data paper",
    tabPanel("Barplot", barpage),
    tabPanel("Rarefaction curves", rarefpage)
    
)

# Define server logic
server <- function(input, output) {
    
    #Barplot panel
    source("panels/panel_server_barplot_dp.R", local = TRUE)
    
    #Rarefaction panel
    source("panels/panel_server_rarefaction_dp.R", local = TRUE)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
