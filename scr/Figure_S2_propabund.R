library(shiny)
library(here)
library(ggplot2)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)
library(mvabund)
library(RColorBrewer)
library(reshape2)
library(purrr)
library(dplyr)
library(plotly)
library(indicspecies)
library(adespatial)
library(fossil)
library(ggdendro)
library(dendextend)



otutab_sfsep0 <- read.table(here("data", "otutab_3.200all_mtax.txt"), header = T, sep = "\t")
otutab_sfsep1 <- cbind.data.frame(otutab_sfsep0[,c(1:8,11,10,9,12:130,133,132,131,134:177)])
names(otutab_sfsep1) <- names(otutab_sfsep0)
otutab_sfsep <- otutab_sfsep1[,-c(141:166)]

pico <- readLines(here("data","pico_descnames.txt"))
three10 <- readLines(here("data","three10_descnames.txt"))
three180 <- readLines(here("data","three180_descnames.txt"))
ten50 <- readLines(here("data","ten50_descnames.txt"))
fifty200 <- readLines(here("data","fifty200_descnames.txt"))
three200 <- readLines(here("data", "three200_descnames.txt"))
three200.2 <- three200[!three200 %in% three180]
net_all <- readLines(here("data","net_all_descnames.txt"))
sfnames <- list("sf0.4.3" = pico, "sf3.10" = three10, "sf10.50" = ten50, "sf50.200" = fifty200, "sf3.180" = three180, "sf3.200" = three200)
sfs <- names(sfnames)

npico <- length(pico)
nthree10 <- length(three10)
nten50 <- length(ten50)
nfifty200 <- length(fifty200)
nthree180 <- length(three180)
nthree200 <- length(three200)

#### Make legend ####
bpcol <- read_delim(here("data","col_figS2.txt"), col_names = T, delim = "\t", comment = "")
bpcol <- bpcol %>% mutate_at(vars(Divisionlong), list( ~factor(., levels = unique(bpcol$Divisionlong), ordered = T)))
bpcolvec <- as.character(bpcol$col)
names(bpcolvec) <- bpcol$Divisionlong

bp_leg_plot <- ggplot(bpcol, aes(x = Divisionlong, fill = Divisionlong))+
  geom_bar()+
  scale_fill_manual(values = bpcolvec)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)),
         fill = guide_legend(nrow = 4, byrow = F))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))
bp_leg_plot
bp_legnd <- get_legend(bp_leg_plot)




taxlevel <- "Divisionlong"



otutab_tax_num <- otutab_sfsep %>% select_if(is.numeric)
####### hei ###
# if (input$propchoice == "selectgroup") {
#   otutab_tax_selectprop <- sweep(otutab_tax_num, 2 , colSums(otutab_tax_num), FUN = "/")
# } else {
  otutab_tax_selectprop <- otutab_tax_num
# }

#####hei###
otutab_tax_notnum <- otutab_sfsep %>% select_if(negate(is.numeric))

otutab_tax <- cbind.data.frame(otutab_tax_notnum, otutab_tax_selectprop)
otutab_tax[is.na(otutab_tax)] <- 0 #NA's because some samples have 0% of a certain class 

#  
otutab_tax$total <- otutab_tax %>% select_if(is.numeric) %>% rowSums()
taxlevtab <- otutab_tax %>% select(get(taxlevel)) %>% droplevels.data.frame()

########

taxgroupspre <-  select(otutab_tax, -total) %>%   group_by_(taxlevel) %>% summarise_if(is.numeric, sum)






limfun <- function(x) {
  ifelse(x>=0.05,1,0)
}

taxgroupspre_mat <- as.matrix(taxgroupspre[,-1, drop = FALSE])
taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun)

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

#if (input$propchoice == "selectgroup") {
taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
#   taxgroups_select2 <- sweep(taxgroups_select2, 2 , colSums(taxgroups_select2), FUN = "/")
# } else {
#   taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
# }
#



taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)

taxgroups_select3t <- t(taxgroups_select3)
taxgroups_select3t2 <- taxgroups_select3t[-1,]
taxgroups_select3t3 <- apply(taxgroups_select3t2, 2, as.numeric)
taxgroups_select3tdf0 <- as.data.frame(taxgroups_select3t3)
names(taxgroups_select3tdf0) <- Taxonomic_group

total04.3 <- taxgroups_select3 %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
total3.10 <- taxgroups_select3 %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
total10.50 <- taxgroups_select3 %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
total50.200<- taxgroups_select3 %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
total3.180 <- taxgroups_select3 %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))


totaltotal04.3 <- otutab_sfsep %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
totaltotal3.10 <- otutab_sfsep %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
totaltotal10.50 <- otutab_sfsep %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
totaltotal50.200<- otutab_sfsep %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
totaltotal3.180 <- otutab_sfsep %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))

taxgroups_select3tdf <- rbind(taxgroups_select3tdf0)

#### 11.07.20 her er jeg ####

taxgroups_select3tdf$sf <- factor(c(rep("sf0.4.3", npico), rep("sf3.10", nthree10),
                                    rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180)), ordered = T,
                                  levels = c("sf0.4.3", "sf3.10", "sf3.180", "sf10.50", "sf50.200"))

taxgroups_select3tdf$stdep <-  factor(c(as.character(StationDepthsfsep_tb)), levels = levels(StationDepthsfsep_tb), ordered = T)
taxgroups_select3tdf$cruise <-  factor(c(cruisesfsep_tb$cruise), ordered = T)


# taxgroups_select3tdf$sf <- factor(c(rep("sf0.4.3", npico), rep("sf3.10", nthree10),
#                      rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180)), ordered = T,
#                      levels = c("sf0.4.3", "sf3.10", "sf3.180", "sf10.50", "sf50.200"))
#
# taxgroups_select3tdf$stdep <-  StationDepthsfsep_tb
# taxgroups_select3tdf$cruise <-  factor(cruisesfsep_tb$cruise, ordered = T)
#


taxgroups_select3tdf_mp12 <- taxgroups_select3tdf %>% filter(cruise %in% c("MP1", "MP2"))

taxgroups_select3tdf_mp345 <- taxgroups_select3tdf %>% filter(cruise %in% c("MP3", "MP4", "MP5"))

taxgroups_select3tdf_total <- taxgroups_select3tdf %>% filter(cruise %in% c("total"))


taxgroups_select3tdf_mp12melt <- melt(taxgroups_select3tdf_mp12)
taxgroups_select3tdf_mp345melt <- melt(taxgroups_select3tdf_mp345)

taxgroups_select3tdf_mp12melt$cruise <- recode(taxgroups_select3tdf_mp12melt$cruise, "MP1" = "Jan", "MP2" = "March") 
taxgroups_select3tdf_mp345melt$cruise <- recode(taxgroups_select3tdf_mp345melt$cruise, "MP3" = "May", "MP4" = "Aug", "MP5" = "Nov")

p12 <- ggplot(taxgroups_select3tdf_mp12melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
  labs(title = "Proportional read abundance")+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ylim(0,1.01)+
  coord_flip()
p12

p345 <- ggplot(taxgroups_select3tdf_mp345melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = bpcolvec)+
  facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ylim(0,1.01)+
  coord_flip()
p345

figureS2 <- ggdraw()+
  draw_plot(p12+theme(legend.position = "none"), x = 0, y = 3/5, width = 0.5, height = 2/5)+
  draw_plot(p345+labs(fill = "Division"), x = 0, y = 0, width = 1, height = 3/5)
figureS2
