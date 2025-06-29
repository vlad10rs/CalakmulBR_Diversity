'Script for abundances graphs and diversity analyses
Paper title: Spatiotemporal diversity of terrestrial mammals in the Maya Forest, Mexico
	Paper authors: Contreras-Moreno, F., et al.
		Script authors: J. Vladimir Rojas Sánchez and J. Juan Flores Martínez'

pacman::p_load(terra, sf, dplyr, mapview, geodata, rgbif, stars, iNEXT, ggplot2, tidyr, ggridges, plotly, ggbump, betapart, boot, ggspatial)
#Abundance graphs
datos <- read.csv("Abundances.csv") ; head(datos)
datos <- datos %>% rename("2016" = X2016, "2017" = X2017, "2018" = X2018, "2019" = X2019, "2020" = X2020, "2021" = X2021, "2022" = X2022)
datos_long <- datos %>% pivot_longer(cols = '2016':'2022', names_to = "Año", values_to = "Abundancia") %>% mutate(Año = as.numeric(Año))

facet_plot <- ggplot(datos_long, aes(x = Año, y = Abundancia, color = Zona, group = interaction(Especie, Zona))) +
	geom_line() + geom_point() + facet_wrap(~ Especie, scales = "free_y") + theme_minimal() + 
	labs(title = "Tendencia de abundancia de especies por año", x = "Año", y = "Abundancia") ; facet_plot

datos_long_ok <- datos_long[which(datos_long$Especie != "Ateles geoffroyi" & datos_long$Especie != "Canis latrans" & datos_long$Especie != "Conepatus semistriatus" &
																		datos_long$Especie != "Herpailurus yagouaroundi" & datos_long$Especie != "Mustela frenata" & datos_long$Especie != "Spilogale gracilis" &
																		datos_long$Especie != "Philander opossum" & datos_long$Especie != "Spilogale yucatanensis" & datos_long$Especie != "Sciurus yucatanensis"),]

datos_long_ok <- rename(datos_long_ok, Zone = Zona)
datos_long_ok <- datos_long_ok %>% mutate(Zone = case_when(Zone == "CBR Inside" ~ "CBR Inside", Zone == "CBR Outside" ~ "CBR Outside"))

scatter_plot <- ggplot(datos_long_ok, aes(x = Año, y = Abundancia, color = Zone)) + geom_point() + geom_smooth(se = FALSE) + 
	facet_wrap(~ Especie, scales = "free_y") + theme_minimal() + labs(x = "Year", y = "Relative abundance index") +
	scale_color_manual(values = c("CBR Inside" = "#3AB7C3", "CBR Outside" = "#F16736")) + 
	theme_minimal() + theme(legend.position = "bottom"); scatter_plot

#Rank abundance graph
datos_rank <- datos_long %>% group_by(Zona, Año) %>% arrange(desc(Abundancia)) %>% mutate(ranking = row_number()) %>% ungroup()

datos_abgroup <- datos_long %>% group_by(Especie, Zona) %>% summarize(mean_abund = mean(Abundancia, na.rm = TRUE)) %>%
	ungroup() %>% mutate(grupo = cut(mean_abund, breaks = quantile(mean_abund, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE),
																	 include.lowest = TRUE, labels = c("Rare", "Common", "Dominant")))

datos_grupos <- merge(datos_long, datos_abgroup[, c("Especie", "Zona", "grupo")], by = c("Especie", "Zona"))

rangab_plot <- ggplot(datos_grupos, aes(x = Año, y = Abundancia, group = Especie, color = grupo)) + geom_line(alpha = 0.7) + geom_point() +
	facet_wrap(~ Zona) +	theme_minimal() + labs(title = "", x = "Year", y = "Abundance", color = "Group") + 
	scale_color_manual(values = c("Rare" = "#F16736", "Common" = "grey", "Dominant" = "#3AB7C3")) ; rangab_plot

#Diversity analyses per Zones
diversityZones <- read.csv("DiversityZones.csv")
diversityZones

divZ <- iNEXT(diversityZones, q = c(0,1,2), datatype = "abundance")
divZ$DataInfo
data.info <- divZ$AsyEst ; data.info
divZ$iNextEst$coverage_based

ggiNEXT(divZ, type = 2, facet.var = "Assemblage", color.var = "Order.q")

divZDF <- data.frame(x = data.info$Assemblage, y = data.info$Observed) ; divZDF
evenness_norm <- function(site) {
	even_n <- divZDF[which(divZDF$x == site),]
	factorEn <- ((even_n[3,2])-1)/((even_n[1,2])-1)
	print(even_n)
	print(factorEn)}
evenness_norm("CBR_Outside") ; evenness_norm("CBR_Inside")

#Diversity analyses per Year
	#Inside the CBR
div_rbc <- read.csv("DiversityRBC.csv") ; div_rb <- iNEXT(div_rbc, q = c(0,1,2), datatype = "abundance")
div_rb$DataInfo
data.info.rb <- div_rb$AsyEst ; data.info.rb ; metrics <- c("Species richness", "Shannon diversity", "Simpson diversity")
dfRBC <- data.frame(Year = data.info.rb$Assemblage, Value = data.info.rb$Observed) 
dfRBC$Metric <- rep(metrics, times = length(unique(dfRBC$Year))) ; dfRBC
dfRBC$Year <- gsub("X", "", dfRBC$Year) ; dfRBC

evenness_norm <- function(site) {
	even_n <- dfRBC[which(dfRBC$Year == site),]
	factorEN <- ((even_n[3,2])-1)/((even_n[1,2])-1)
	print(even_n)
	print(factorEN)}
evenness_norm("2016") ;  evenness_norm("2017") ; evenness_norm("2018"); evenness_norm("2019") ; evenness_norm("2020"); evenness_norm("2021") ; evenness_norm("2022")

	#Outside the CBR
div_out <- read.csv("DiversityOutside.csv") ; div_out <- iNEXT(div_out, q = c(0,1,2), datatype = "abundance")
div_out$DataInfo
data.info.out <- div_out$AsyEst ; data.info.out ; metrics <- c("Species richness", "Shannon diversity", "Simpson diversity")
dfOut <- data.frame(Year = data.info.out$Assemblage, Value = data.info.out$Observed) 
dfOut$Metric <- rep(metrics, times = length(unique(dfRBC$Year))) ; dfOut
dfOut$Year <- gsub("X", "", dfOut$Year) ; dfOut

evenness_norm <- function(site) {
	even_n <- dfOut[which(dfOut$Year == site),]
	factorEN <- ((even_n[3,2])-1)/((even_n[1,2])-1)
	print(even_n)
	print(factorEN)}
evenness_norm("2016") ;  evenness_norm("2017") ; evenness_norm("2018"); evenness_norm("2019") ; evenness_norm("2020"); evenness_norm("2021") ; evenness_norm("2022")

df_rbc <- data.frame(Evenness = c(0.2571, 0.4054, 0.2767, 0.3104, 0.3718, 0.3571, 0.3408), Year = c(2016:2022))
df_out <- data.frame(Evenness = c(0.5691, 0.4080, 0.5784, 0.5581, 0.3815, 0.5002, 0.4358), Year = c(2016:2022))

dual_line_plot <- ggplot(df_combined, aes(x = Year_numeric, y = Evenness, color = Source, group = Source)) + geom_point(size = 2) +
	geom_smooth(method = loess, se = FALSE, size = 1, linetype = "solid") +
	theme_minimal() + labs(x = "Year", y = "Simpson Evenness value", color = "Zone") +
	scale_color_manual(values = c("CBR inside" = "#3AB7C3", "CBR outside" = "#F16736")) + theme_minimal() + theme(legend.position = "bottom") ; dual_line_plot
