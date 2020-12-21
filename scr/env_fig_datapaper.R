dt <- tibble("month" = c(rep("Jan",4), rep("Mar",4)), "st" =c("s1", "s1", "s2","s2","s3","s3","s4","s4"), "depth" = c(1,9,1,10,1,12,1,11), nut = c(0.1,0.2,0.5,.6, .4, .5, .2, .1))
meta_table <- read.table(here("data", "meta_cleanMP.txt"), header = T, sep = "\t", fill = NA)
profiletab <- read.table(here("data", "MP_CB_profiles_20200206.txt"), header = T, sep = "\t", fill = NA, na.strings = c("", " ", "NA"))
head(profiletab)

profiletab2 <- profiletab[rowSums(is.na(profiletab)) != ncol(profiletab),]
head(profiletab2)
names(profiletab2)

meta_table$Cruise <- recode(meta_table$Cruise, "MP1" = "Jan", "MP2" = "March", "MP3" = "May", "MP4" = "Aug", "MP5" = "Nov") 
profiletab2$Cruise <- recode(profiletab2$Cruise, "MP1" = "Jan", "MP2" = "March", "MP3" = "May", "MP4" = "Aug", "MP5" = "Nov") 

#### make linetype for each station ####
(text=expression(paste("Chl ", italic(a), " [", mu,"g",L^-1,"]", sep="")),  side=2, line=2)
#mtext(text=expression(paste("Temperature [", degree,"C]      PAR [mol ", m^-2, "da",y^-1,"]       Salinity [PSU]", sep="")), side=2, line=2.5, cex=0.8)
mtext(text=expression(paste("[", degree,"C] [mol ", m^-2, "da",y^-1,"] [PSU]", sep="")), side=2, line=2.5, cex=0.8)

axte <- 0.6
axti <- 0.6       
lts <- tibble("Station" = levels(meta_table$Station), "lityp" = c(2,1,1,2,3,4,5,1,2,3,1,2,3,1,2,3))


meta_phys <- profiletab2 %>% select(Cruise, Station, DEPTH_M, CTD.S, CTD.T, sigma)  %>% mutate(DEPTH_M = as.numeric(as.character(DEPTH_M))) %>%  filter(Station %in% unique(lts$Station) & DEPTH_M<= 1000)
meta_phys <- left_join(meta_phys, lts, by ="Station") %>% mutate(lityp = as.factor(lityp))
meta_phys_melt <- melt(meta_phys, measure.vars = c("CTD.S", "CTD.T", "sigma")) #"variable" is now the type of parameter
meta_phys_melt

meta_phys_points <- meta_table %>% select(Cruise, Station, DEPTH_M, CTD.S, CTD.T, sigma)  %>% mutate(DEPTH_M = as.numeric(as.character(DEPTH_M))) %>%  filter(Station %in% unique(lts$Station) & DEPTH_M<= 1000)
meta_phys_points <- left_join(meta_phys_points, lts, by ="Station") %>% mutate(lityp = as.factor(lityp))
meta_phys_points_melt <- melt(meta_phys_points, measure.vars = c("CTD.S", "CTD.T", "sigma")) #"variable" is now the type of parameter
meta_phys_points_melt

sal <- ggplot(data = meta_phys_melt %>% filter(variable == "CTD.S")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_phys_points_melt%>% filter(variable == "CTD.S"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Salinity [PSU]")+
  coord_flip()
sal

temp <- ggplot(data = meta_phys_melt %>% filter(variable == "CTD.T")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_phys_points_melt%>% filter(variable == "CTD.T"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Temperature [", degree,"C]")))+
  coord_flip()
temp

sigma <- ggplot(data = meta_phys_melt %>% filter(variable == "sigma")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_phys_points_melt%>% filter(variable == "sigma"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P5, P4, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Density ", sigma)))+
  coord_flip()
sigma

meta_chem <- profiletab2 %>% select(Cruise, Station, DEPTH_M, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>% mutate(DEPTH_M = as.numeric(as.character(DEPTH_M))) %>%  filter(Station %in% unique(lts$Station) & DEPTH_M<= 1000)
meta_chem <- left_join(meta_chem, lts, by ="Station") %>% mutate(lityp = as.factor(lityp))
meta_chem_melt <- melt(meta_chem, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter
meta_chem_melt

meta_chem_points <- meta_table %>% select(Cruise, Station, DEPTH_M, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)
meta_chem_points <- left_join(meta_chem_points, lts, by ="Station") %>% mutate(lityp = as.factor(lityp))
meta_chem_points_melt <- melt(meta_chem_points, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a"))

inorgN <- ggplot(data = meta_chem_melt %>% filter(variable == "Total_inorgN_uM")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_chem_points_melt%>% filter(variable == "Total_inorgN_uM"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Total inorganic N", " [", mu,"M]")))+
  coord_flip()
inorgN

phos <- ggplot(data = meta_chem_melt %>% filter(variable == "PO4_uM")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_chem_points_melt%>% filter(variable == "PO4_uM"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("PO"[4]^"3-", " [", mu,"M]")))+
  coord_flip()
phos

sil <- ggplot(data = meta_chem_melt %>% filter(variable == "SiOH4_uM")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_chem_points_melt%>% filter(variable == "SiOH4_uM"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("SiOH"[4], " [", mu,"M]")))+
  coord_flip()
sil

chl <- ggplot(data = meta_chem_melt %>% filter(variable == "Chl_a")) + aes(x = desc(sqrt(DEPTH_M)), y = value, lty = lityp)+
  geom_point(data = meta_chem_points_melt%>% filter(variable == "Chl_a"), aes(x = desc(sqrt(DEPTH_M)), y = value, pch = lityp))+
  geom_line()+
  facet_grid(.~ Cruise)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)),
        axis.title = element_text(size = rel(axti)),
        axis.title.y = element_blank())+
  xlab("Depth")+
  ylab(expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))+
  coord_flip()
chl

envfig <- ggarrange(sal+theme(legend.title =element_text(size = 8),
                   legend.text = element_text(size = 8), axis.title.y = element_blank()), temp+theme(axis.title.y = element_blank()), sigma+theme(axis.title.y = element_blank()), 
                   inorgN+theme(axis.title.y = element_blank()), phos+theme(axis.title.y = element_blank()), sil+theme(axis.title.y = element_blank()), chl+theme(axis.title.y = element_blank()), nrow = 7, common.legend = TRUE, legend = "right")
  
annotate_figure(envfig, left = text_grob("Depth [m]", rot = 90))

