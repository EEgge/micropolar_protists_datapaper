output$salplotly <- renderPlotly({
  env_data_figure()$salplotly
})

output$tempplotly <- renderPlotly({
  env_data_figure()$tempplotly
})

output$sigmaplotly <- renderPlotly({
  env_data_figure()$sigmaplotly
})

output$inorgNplotly <- renderPlotly({
  env_data_figure()$inorgNplotly
})

output$phosplotly <- renderPlotly({
  env_data_figure()$phosplotly
})

output$silplotly <- renderPlotly({
  env_data_figure()$silplotly
})

output$chlplotly <- renderPlotly({
  env_data_figure()$chlplotly
})

env_data_figure <- reactive({

axte <- 0.6
axti <- 0.8       
lts <- tibble("station" = levels(as.factor(env_table$station)), "lityp" = c(1,2,1,2,3,4,5,1,2,3,1,2,3,1,2,3))

env_profiles <- env_profiles %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_phys <- env_profiles %>% select(month, station, depth_m, CTD.S, CTD.T, sigma)  %>%  filter(depth_m <= 1000)
env_phys <- left_join(env_phys, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_phys_melt <- reshape2::melt(env_phys, measure.vars = c("CTD.S", "CTD.T", "sigma")) #"variable" is now the type of parameter
env_phys_melt


env_table <- env_table %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_phys_points <- env_table %>% select(month, station, depth_m, CTD.S, CTD.T, sigma) %>%  filter(depth_m <= 1000)
env_phys_points <- left_join(env_phys_points, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_phys_points_melt <-reshape2::melt(env_phys_points, measure.vars = c("CTD.S", "CTD.T", "sigma")) #"variable" is now the type of parameter
env_phys_points_melt

salplot <- ggplot(data = env_phys_melt %>% filter(variable == "CTD.S"), aes(x = desc(sqrt(depth_m)), y = value, lty = lityp))+
  geom_line()+
  geom_point(data = env_phys_points_melt%>% filter(variable == "CTD.S"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>PSU: %s", station, depth_m, round(value,2))))+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  scale_shape_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Salinity [PSU]")+
  coord_flip()
salplot
salplotly <- ggplotly(salplot, tooltip = "text")

tempplot <- ggplot(data = env_phys_melt %>% filter(variable == "CTD.T")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_phys_points_melt%>% filter(variable == "CTD.T"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>DegC: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Temperature [deg C]")+
  coord_flip()
tempplot
tempplotly <- ggplotly(tempplot, tooltip = "text")
tempplotly

sigmaplot <- ggplot(data = env_phys_melt %>% filter(variable == "sigma")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_phys_points_melt%>% filter(variable == "sigma"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>Sigma: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P5, P4, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Density (sigma)")+
  coord_flip()
sigmaplot
sigmaplotly <- ggplotly(sigmaplot, tooltip = "text")

env_chem <- env_profiles %>% select(month, station, depth_m, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a) %>%  filter(depth_m<= 1000)
env_chem <- left_join(env_chem, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_chem_melt <- reshape2::melt(env_chem, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter
env_chem_melt

env_chem_points <- env_table %>% select(month, station, depth_m, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)
env_chem_points <- left_join(env_chem_points, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_chem_points_melt <- reshape2::melt(env_chem_points, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a"))

inorgNplot <- ggplot(data = env_chem_melt %>% filter(variable == "Total_inorgN_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_chem_points_melt%>% filter(variable == "Total_inorgN_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("Total inorganic N [microM]")+
  coord_flip()
inorgNplot
inorgNplotly <- ggplotly(inorgNplot, tooltip = "text")

phosplot <- ggplot(data = env_chem_melt %>% filter(variable == "PO4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_chem_points_melt%>% filter(variable == "PO4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("PO43-, [microM]")+
  coord_flip()
phosplot
phosplotly <- ggplotly(phosplot, tooltip = "text")

silplot <- ggplot(data = env_chem_melt %>% filter(variable == "SiOH4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_chem_points_melt%>% filter(variable == "SiOH4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab("SiOH4, [microM]")+
  coord_flip()
silplot
silplotly <- ggplotly(silplot, tooltip = "text")

chlplot <- ggplot(data = env_chem_melt %>% filter(variable == "Chl_a")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
  geom_point(data = env_chem_points_melt%>% filter(variable == "Chl_a"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microgL-1: %s", station, depth_m, round(value,2))))+
  geom_line()+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
  theme(axis.text = element_text(size = rel(axte)),
        axis.title = element_text(size = rel(axti)),
        axis.title.y = element_blank())+
  xlab("Depth")+
  ylab("Chl a, [ microgL^-1]")+
  coord_flip()
chlplot
chlplotly <- ggplotly(chlplot, tooltip = "text")

envfig <- ggarrange(salplot+theme(legend.title =element_text(size = 8),
                              legend.text = element_text(size = 8), axis.title.y = element_blank()), tempplot+theme(axis.title.y = element_blank()), sigmaplot+theme(axis.title.y = element_blank()), 
                    inorgNplot+theme(axis.title.y = element_blank()), phosplot+theme(axis.title.y = element_blank()), silplot+theme(axis.title.y = element_blank()), chlplot+theme(axis.title.y = element_blank()), nrow = 7, common.legend = TRUE, legend = "right")

annotate_figure(envfig, left = text_grob("Depth [m]", rot = 90, size = 8))

list(salplotly = salplotly, tempplotly = tempplotly, sigmaplotly = sigmaplotly, inorgNplotly = inorgNplotly, phosplotly = phosplotly, silplotly = silplotly, chlplotly = chlplotly)

})
