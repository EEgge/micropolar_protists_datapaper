#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Interactive figures from MicroPolar protist metabarcoding data paper"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("prop_lim", "Taxonomic groups present with < x% in all samples will be grouped together as 'Other'",
                         value = 5, min = 0, max = 100),
            radioButtons("which_tab", "ASV-table with separate or merged replicates:", 
                         choices = c("Merged" = "merged", "Separate" = "separate"), selected = "merged"),
            checkboxGroupInput("taxo_group2", label = "Division", 
                               choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                           "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Pseudofungi",  "Radiolaria",  "All"), 
                               selected = "All"),
            radioButtons("taxlevel_barplot", label = "Taxonomic level",
                         choices = c("Division" = "Divisionlong", "Class"="Class", "Order"="Order", "Family"="Family", "Genus"="Genus", "Species"="Species")),
            radioButtons("propchoice", label = "Proportion of selected group or all", choices = c("Selected group" = "selectgroup","All" = "All"))
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           textOutput("jaja")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    limfun <- reactive({
        function(x) {
            ifelse(x>=input$prop_lim/100,1,0)
        } 
    })
    
    mp_plot <- reactive({ 
        taxlevel <- input$taxlevel_barplot
        otutab_sfsep <- otutab_sfsep %>% select(-all_of(three200.2)) #Remove 'fake' size fraction 3-200
        
        if (input$taxo_group2 != "All") {
            otutab_tax0 <- otutab_sfsep %>% filter(.data[[input$taxlevel]] %in% .env[[input$taxo_group2]]) #Create new OTU table with only the selected division
        } else {
            otutab_tax0 <- otutab_sfsep
        }
        
        otutab_tax_num <- otutab_tax0 %>% select_if(is.numeric)
        
        #### If selected, create new proportions of selected groups ####
        otutab_tax_num <- otutab_sfsep %>% select_if(is.numeric)
        if (input$propchoice == "selectgroup") {
          otutab_tax_selectprop <- sweep(otutab_tax_num, 2 , colSums(otutab_tax_num), FUN = "/")
        } else {
        otutab_tax_selectprop <- otutab_tax_num
         }
        
        
        otutab_tax_notnum <- otutab_sfsep %>% select_if(negate(is.numeric))
        otutab_tax <- cbind.data.frame(otutab_tax_notnum, otutab_tax_selectprop)
        otutab_tax[is.na(otutab_tax)] <- 0 #NA's because some samples have 0% of a certain class 
        
        #  
        otutab_tax$total <- otutab_tax %>% select_if(is.numeric) %>% rowSums()
        taxlevtab <- otutab_tax %>% select(get(taxlevel)) %>% droplevels.data.frame()
        
        ########
        
        taxgroupspre <-  select(otutab_tax, -total) %>%   group_by_(taxlevel) %>% summarise_if(is.numeric, sum)
        
        
        
        
        
        
        
        
        taxgroupspre_mat <- as.matrix(taxgroupspre[,-1, drop = FALSE])
        taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun())
        
        if (is.vector(taxgroupspre_bin)) {
            taxgroupspre_bin2 <- as.data.frame(as.list(taxgroupspre_bin))
        } else {
            taxgroupspre_bin2 <- as.data.frame(taxgroupspre_bin)}
        
        
        
        taxgroupspre_bin_sums <- taxgroupspre_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
        taxgroupspre_bin_sums$taxgroups <- as.vector(taxgroupspre[[taxlevel]])
        taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% select(taxgroups)
        
        
        
        
        taxgroups_select <- taxgroupspre %>% filter(get(taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
        
        taxgroups_other <- taxgroupspre %>% filter(!get(taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
        
        
        
        if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
            Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
            othersum <- colSums(taxgroups_other[,-1])
        } else {
            Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
            othersum <- NULL
        }
        #sjekk_prop <- head(taxgroupspre_mat)[1,] ####14.03.2020
        
        if (input$propchoice == "selectgroup") {
        taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
           taxgroups_select2 <- sweep(taxgroups_select2, 2 , colSums(taxgroups_select2), FUN = "/")
         } else {
           taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
         }
        
        
        
        
        taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)
        
        taxgroups_select3t <- t(taxgroups_select3)
        taxgroups_select3t2 <- taxgroups_select3t[-1,]
        taxgroups_select3t3 <- apply(taxgroups_select3t2, 2, as.numeric)
        taxgroups_select3tdf0 <- as.data.frame(taxgroups_select3t3)
        names(taxgroups_select3tdf0) <- Taxonomic_group
        
        
    })
    
    output$jaja <- renderText({
        limfun()(c(.06,.02))
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
