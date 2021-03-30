
env_table <- read_delim(here("data", "env_data_depths.txt"), delim = "\t")
env_profiles <- read_delim(here("data", "env_data_profiles.txt"), delim = "\t")

axte <- 0.8
axti <- 1
annsize <- 8
lts <- tibble("station" = levels(as.factor(env_table$station)), "lityp" = c(2,1,1,2,3,4,23,1,2,3,1,2,3,1,2,3))

env_profiles <- env_profiles %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_profiles <- env_profiles %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>%  filter(depth_m <= 1000)
env_profiles <- left_join(env_profiles, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_profiles_melt <- reshape2::melt(env_profiles, measure.vars = c("CTD.S", "CTD.T", "sigma","Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter

env_points <- env_table %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_points <- env_points %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>% mutate(depth_m = as.numeric(as.character(depth_m))) %>%  filter(station %in% unique(lts$station) & depth_m<= 1000)
env_points <- left_join(env_points, lts, by ="station") %>% mutate(lityp = as.character(lityp))
env_points_melt <- reshape2::melt(env_points, measure.vars = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter

envplot <- function(envvar_profile, envvar_point, ylab) {
  ggplot(data = envvar_profile) + 
    aes(x = desc(sqrt(depth_m)), y = value, lty = lityp, col =lityp, group = lityp)+
    geom_point(data = envvar_point, aes(x = desc(sqrt(depth_m)), y = value, shape = lityp, col =lityp))+
    geom_line()+
    facet_grid(.~ month)+
    theme_classic()+
    theme(strip.background = element_rect(colour="white", fill="white"))+
    scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
    scale_linetype_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
    scale_color_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
    scale_shape_manual(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"), values = c(19,17,3,4,7))+
    theme(axis.text = element_text(size = rel(axte)))+
    theme(axis.title = element_text(size = rel(axti)),
          strip.text = element_text(size = 11))+
    xlab("Depth")+
    ylab(ylab)+
    coord_flip()+
    theme(legend.position = "right",
          legend.text = element_text(size = 13))+
    #guides(shape = guide_legend(override.aes = list(size = 3)))+
    theme(legend.key.size = unit(2,"line"))
}


sal <- envplot(env_profiles_melt %>% filter(variable == "CTD.S"), env_points_melt%>% filter(variable == "CTD.S"), "Salinity [PSU]")
sal
temp <- envplot(env_profiles_melt %>% filter(variable == "CTD.T"), env_points_melt%>% filter(variable == "CTD.T"), expression(paste("Temperature [", degree,"C]")))
sigma <- envplot(env_profiles_melt %>% filter(variable == "sigma"), env_points_melt%>% filter(variable == "sigma"), expression(paste("Density ", sigma)))
inorgN <- envplot(env_profiles_melt %>% filter(variable == "Total_inorgN_uM"), env_points_melt%>% filter(variable == "Total_inorgN_uM"), expression(paste("Total inorganic N", " [", mu,"M]")))
phos <- envplot(env_profiles_melt %>% filter(variable == "PO4_uM"), env_points_melt%>% filter(variable == "PO4_uM"), expression(paste("PO"[4]^"3-", " [", mu,"M]")))
sil <- envplot(env_profiles_melt %>% filter(variable == "SiOH4_uM"), env_points_melt%>% filter(variable == "SiOH4_uM"), expression(paste("SiOH"[4], " [", mu,"M]")))
chl <- envplot(env_profiles_melt %>% filter(variable == "Chl_a"), env_points_melt%>% filter(variable == "Chl_a"), expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))

envfig <- ggarrange(sal+theme(legend.title =element_text(size = 10),
                   legend.text = element_text(size =9), axis.title.y = element_blank()), temp+theme(axis.title.y = element_blank()), sigma+theme(axis.title.y = element_blank()), 
                   inorgN+theme(axis.title.y = element_blank()), phos+theme(axis.title.y = element_blank()), sil+theme(axis.title.y = element_blank()), chl+theme(axis.title.y = element_blank()), nrow = 7, common.legend = TRUE, legend = "right",
                   labels = LETTERS[1:7])
envfig

  
