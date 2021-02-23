dt <- tibble("month" = c(rep("Jan",4), rep("Mar",4)), "st" =c("s1", "s1", "s2","s2","s3","s3","s4","s4"), "depth" = c(1,9,1,10,1,12,1,11), nut = c(0.1,0.2,0.5,.6, .4, .5, .2, .1))

env_table <- read_delim(here("data", "env_data_depths.txt"), delim = "\t")
env_profiles <- read_delim(here("data", "env_data_profiles.txt"), delim = "\t")

axte <- 0.7
axti <- 0.7
annsize <- 8
lts <- tibble("station" = levels(as.factor(env_table$station)), "lityp" = c(2,1,1,2,3,4,5,1,2,3,1,2,3,1,2,3))

env_profiles <- env_profiles %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_profiles <- env_profiles %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>%  filter(depth_m <= 1000)
env_profiles <- left_join(env_profiles, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_profiles_melt <- reshape2::melt(env_profiles, measure.vars = c("CTD.S", "CTD.T", "sigma","Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter
env_profiles_melt

env_points <- env_table %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_points <- env_points %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>% mutate(depth_m = as.numeric(as.character(depth_m))) %>%  filter(station %in% unique(lts$station) & depth_m<= 1000)
env_points <- left_join(env_points, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_points_melt <- reshape2::melt(env_points, measure.vars = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter
env_points_melt

sal <- ggplot(data = env_profiles_melt %>% filter(variable == "CTD.S")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt%>% filter(variable == "CTD.S"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  scale_shape_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  scale_color_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Salinity [PSU]")+
  coord_flip()+
  theme(legend.position = "top")+
  #guides(shape = guide_legend(override.aes = list(size = 1)))
theme(legend.key.size = unit(2,"line"))
sal
sal2 <- annotate_figure(sal, right = text_grob("Salinity", rot = 0, size = annsize))

temp <- ggplot(data = env_profiles_melt %>% filter(variable == "CTD.T")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt%>% filter(variable == "CTD.T"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Temperature [", degree,"C]")))+
  coord_flip()+
  theme(legend.position = "top")
temp
temp2 <- annotate_figure(temp, right = text_grob("Temperature", rot = 0, size = annsize))

sigma <- ggplot(data = env_profiles_melt %>% filter(variable == "sigma")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt%>% filter(variable == "sigma"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P5, P4, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Density ", sigma)))+
  coord_flip()+
  theme(legend.position = "top")
sigma
sigma2 <- annotate_figure(sigma, right = text_grob("Density", rot = 0, size = annsize))

inorgN <- ggplot(data = env_profiles_melt %>% filter(variable == "Total_inorgN_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt%>% filter(variable == "Total_inorgN_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("Total inorganic N", " [", mu,"M]")))+
  coord_flip()+
  theme(legend.position = "top")
inorgN
inorgN2 <- annotate_figure(inorgN, right = text_grob("Total inorganic N", rot = 0, size = annsize))


phos <- ggplot(data = env_profiles_melt %>% filter(variable == "PO4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt %>% filter(variable == "PO4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("PO"[4]^"3-", " [", mu,"M]")))+
  coord_flip()+
  theme(legend.position = "top")
phos
phos2 <- annotate_figure(phos, right = text_grob(expression(paste("PO"[4]^"3-")), rot = 0, size = annsize))


sil <- ggplot(data = env_profiles_melt %>% filter(variable == "SiOH4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt %>% filter(variable == "SiOH4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(expression(paste("SiOH"[4], " [", mu,"M]")))+
  coord_flip()+
  theme(legend.position = "top")
sil
sil2 <- annotate_figure(sil, right = text_grob(expression(paste("SiOH"[4])), rot = 0, size = annsize))


chl <- ggplot(data = env_profiles_melt %>% filter(variable == "Chl_a")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp)+
  geom_point(data = env_points_melt%>% filter(variable == "Chl_a"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, col =lityp))+
  geom_line()+
  facet_grid(.~ month)+
  theme_classic()+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)),
        axis.title = element_text(size = rel(axti)),
        axis.title.y = element_blank())+
  xlab("Depth")+
  ylab(expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))+
  coord_flip()+
  theme(legend.position = "top")
chl
chl2 <- annotate_figure(chl, right = text_grob(expression(paste("Chl ",italic(a))), rot = 0, size = annsize))


envfig <- ggarrange(sal+theme(legend.title =element_text(size = 9),
                   legend.text = element_text(size =9), axis.title.y = element_blank()), temp+theme(axis.title.y = element_blank()), sigma+theme(axis.title.y = element_blank()), 
                   inorgN+theme(axis.title.y = element_blank()), phos+theme(axis.title.y = element_blank()), sil+theme(axis.title.y = element_blank()), chl+theme(axis.title.y = element_blank()), nrow = 7, common.legend = TRUE, legend = "top",
                   labels = LETTERS[1:7])
envfig
{phantom("A")~"-"}

  
annotate_figure(envfig, left = text_grob("Depth [m]", rot = 90), right = text_grob(expression(paste('Salinity\n\n\n\n\nTemperature\n\n\n\n\nDensity\n\n\n\n\nTotal inorganic N\n\n\n\n\n'~PO[4]^{"3"~"-"} ~ '\nhm\n\n\n\n\n'~SiOH[4]~'hm \n\n\n\n\nChl'~italic(a))), rot = 0, size = 10))

annotate_figure(envfig, left = text_grob("Depth [m]", rot = 90), right = text_grob(expression(paste('Salinity', 'Temperature', 'Density', 'Total inorganic N', PO[4]^{"3"~"-"}, SiOH[4], sep = "\n")), rot = 0, size = 10))

side_text <- as.list(c("Salinity", "Temperature", "Density", "Total inorganic N"), expression(paste(PO[4]^{"3"~"-"})))

annotate_figure(envfig, right = text_grob(side_text))

side_text2 <- paste("Salinity", "Temperature", expression(paste(SiOH[4])), sep = "\n")

