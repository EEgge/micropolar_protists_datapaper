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
  
  waiter <- waiter::Waiter$new()
  waiter$show()
  on.exit(waiter$hide())
  
  Sys.sleep(sample(5, 1))
  runif(1)

axte <- 0.6
axti <- 0.8       
lts <- tibble("station" = levels(as.factor(env_depths$station)), "lityp" = c(1,2,1,2,3,4,5,1,2,3,1,2,3,1,2,3))

env_profiles <- env_profiles %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_profiles <- env_profiles %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)  %>%  filter(depth_m <= 1000)
env_profiles <- left_join(env_profiles, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_profiles_melt <- reshape2::melt(env_profiles, measure.vars = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter

env_depths <- env_depths %>% mutate(month = factor(month, levels = c("Jan", "Mar", "May", "Aug", "Nov", ordered =T)))
env_depths <- env_depths %>% select(month, station, depth_m, CTD.S, CTD.T, sigma, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a) %>%  filter(depth_m <= 1000)
env_depths <- left_join(env_depths, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
env_depths_melt <-reshape2::melt(env_depths, measure.vars = c("CTD.S", "CTD.T", "sigma", "Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a"))

envplot <- function(envvar_profile, envvar_depths, ylab) {
  envplot <- ggplot(data = envvar_profile, aes(x = desc(sqrt(depth_m)), y = value, lty = lityp))+
  geom_line()+
  geom_point(data = envvar_depths, aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf(paste0("Station: %s<br>Depth: %s<br>", ylab, ": %s"), station, depth_m, round(value,2))))+
  facet_grid(.~ month)+
  scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
  scale_linetype_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  scale_shape_discrete(name = "Station", labels = c("B08, M02, P01, P05, N02", "B16, M03, P03, P06, N03", "M04, P04, P07, N04", "M05", "M06"))+
  theme(axis.text = element_text(size = rel(axte)))+
  theme(axis.title = element_text(size = rel(axti)))+
  xlab("Depth")+
  ylab(ylab)+
  coord_flip()

  envplotly <- ggplotly(envplot, tooltip = "text")
  return(envplotly)
}

salplotly <- envplot(env_profiles_melt %>% filter(variable == "CTD.S"), env_depths_melt%>% filter(variable == "CTD.S"), "Salinity [PSU]")
tempplotly <- envplot(env_profiles_melt %>% filter(variable == "CTD.T"), env_depths_melt%>% filter(variable == "CTD.T"), expression(paste("Temperature [", degree,"C]")))
sigmaplotly <- envplot(env_profiles_melt %>% filter(variable == "sigma"), env_depths_melt%>% filter(variable == "sigma"), expression(paste("Density ", sigma)))
inorgNplotly <- envplot(env_profiles_melt %>% filter(variable == "Total_inorgN_uM"), env_depths_melt%>% filter(variable == "Total_inorgN_uM"), expression(paste("Total inorganic N", " [", mu,"M]")))
phosplotly <- envplot(env_profiles_melt %>% filter(variable == "PO4_uM"), env_depths_melt%>% filter(variable == "PO4_uM"), expression(paste("PO"[4]^"3-", " [", mu,"M]")))
silplotly <- envplot(env_profiles_melt %>% filter(variable == "SiOH4_uM"), env_depths_melt%>% filter(variable == "SiOH4_uM"), expression(paste("SiOH"[4], " [", mu,"M]")))
chlplotly <- envplot(env_profiles_melt %>% filter(variable == "Chl_a"), env_depths_melt%>% filter(variable == "Chl_a"), expression(paste("Chl ",italic(a), " [", mu,"g",L^-1,"]")))


# tempplot <- ggplot(data = env_phys_melt %>% filter(variable == "CTD.T")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_phys_points_melt%>% filter(variable == "CTD.T"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>DegC: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)))+
#   theme(axis.title = element_text(size = rel(axti)))+
#   xlab("Depth")+
#   ylab("Temperature [deg C]")+
#   coord_flip()
# tempplot
# tempplotly <- ggplotly(tempplot, tooltip = "text")
# tempplotly
# 
# sigmaplot <- ggplot(data = env_phys_melt %>% filter(variable == "sigma")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_phys_points_melt%>% filter(variable == "sigma"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>Sigma: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P5, P4, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)))+
#   theme(axis.title = element_text(size = rel(axti)))+
#   xlab("Depth")+
#   ylab("Density (sigma)")+
#   coord_flip()
# sigmaplot
# sigmaplotly <- ggplotly(sigmaplot, tooltip = "text")
# 
# env_chem <- env_profiles %>% select(month, station, depth_m, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a) %>%  filter(depth_m<= 1000)
# env_chem <- left_join(env_chem, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
# env_chem_melt <- reshape2::melt(env_chem, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a")) #"variable" is now the type of parameter
# env_chem_melt
# 
# env_chem_points <- env_table %>% select(month, station, depth_m, Total_inorgN_uM, PO4_uM, SiOH4_uM, Chl_a)
# env_chem_points <- left_join(env_chem_points, lts, by ="station") %>% mutate(lityp = as.factor(lityp))
# env_chem_points_melt <- reshape2::melt(env_chem_points, measure.vars = c("Total_inorgN_uM", "PO4_uM", "SiOH4_uM", "Chl_a"))
# 
# inorgNplot <- ggplot(data = env_chem_melt %>% filter(variable == "Total_inorgN_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_chem_points_melt%>% filter(variable == "Total_inorgN_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)))+
#   theme(axis.title = element_text(size = rel(axti)))+
#   xlab("Depth")+
#   ylab("Total inorganic N [microM]")+
#   coord_flip()
# inorgNplot
# inorgNplotly <- ggplotly(inorgNplot, tooltip = "text")
# 
# phosplot <- ggplot(data = env_chem_melt %>% filter(variable == "PO4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_chem_points_melt%>% filter(variable == "PO4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)))+
#   theme(axis.title = element_text(size = rel(axti)))+
#   xlab("Depth")+
#   ylab("PO43-, [microM]")+
#   coord_flip()
# phosplot
# phosplotly <- ggplotly(phosplot, tooltip = "text")
# 
# silplot <- ggplot(data = env_chem_melt %>% filter(variable == "SiOH4_uM")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_chem_points_melt%>% filter(variable == "SiOH4_uM"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microM: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)))+
#   theme(axis.title = element_text(size = rel(axti)))+
#   xlab("Depth")+
#   ylab("SiOH4, [microM]")+
#   coord_flip()
# silplot
# silplotly <- ggplotly(silplot, tooltip = "text")
# 
# chlplot <- ggplot(data = env_chem_melt %>% filter(variable == "Chl_a")) + aes(x = desc(sqrt(depth_m)), y = value, lty = lityp)+
#   geom_point(data = env_chem_points_melt%>% filter(variable == "Chl_a"), aes(x = desc(sqrt(depth_m)), y = value, pch = lityp, text = sprintf("Station: %s<br>Depth: %s<br>microgL-1: %s", station, depth_m, round(value,2))))+
#   geom_line()+
#   facet_grid(.~ month)+
#   scale_x_continuous(labels = c("0", "20", "500", "1000"), breaks = c(0, -sqrt(20), -sqrt(500), -sqrt(1000)))+
#   scale_linetype_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   scale_shape_discrete(name = "Station", labels = c("B8, M2, P1, P5, N2", "B16, M3, P3, P6, N3", "M4, P4, P7, N4", "M5", "M6"))+
#   theme(axis.text = element_text(size = rel(axte)),
#         axis.title = element_text(size = rel(axti)),
#         axis.title.y = element_blank())+
#   xlab("Depth")+
#   ylab("Chl a, [ microgL^-1]")+
#   coord_flip()
# chlplot
# chlplotly <- ggplotly(chlplot, tooltip = "text")

list(salplotly = salplotly, tempplotly = tempplotly, sigmaplotly = sigmaplotly, inorgNplotly = inorgNplotly, phosplotly = phosplotly, silplotly = silplotly, chlplotly = chlplotly)

})
