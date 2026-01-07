'Script for abundances graphs and diversity analyses
Paper title: Spatiotemporal diversity of terrestrial mammals in the Maya Forest, Mexico
Paper authors: Contreras-Moreno, F., et al.
Script authors: J. Vladimir Rojas Sánchez and J. Juan Flores Martínez'

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, iNEXT, rlang)

#Abundance analysis
datos <- read.csv("Abundances.csv")

#Renaming columns
datos_long <- datos %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "Year", 
               values_to = "Abundance",
               names_prefix = "X") %>%
  mutate(Year = as.numeric(Year))

#Exclude species only present inside the CBR from the plot
exclude_species <- c("Ateles geoffroyi", "Canis latrans", "Conepatus semistriatus",
    "Herpailurus yagouaroundi", "Mustela frenata", "Spilogale gracilis",
    "Philander opossum", "Spilogale yucatanensis", "Sciurus yucatanensis")

datos_long_ok <- datos_long %>%
  filter(!Especie %in% exclude_species) %>%
  rename(Zone = Zona)

scatter_plot <- ggplot(datos_long_ok, aes(x = Year, y = Abundance, color = Zone)) +
  geom_point() +
  geom_smooth(se = FALSE) + 
  facet_wrap(~ Especie, scales = "free_y") +
  scale_color_manual(values = c("CBR Inside" = "#3AB7C3", "CBR Outside" = "#F16736")) + 
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "Relative abundance index") ; scatter_plot

#Diversity analyses
div_zones <- read.csv("div_zones.csv")
div_years <- read.csv("div_years.csv")

#Select the wished scale of analysis
#Options: "zones", "years_inside", "years_outside"
diversity_analysis <- function(scale) {
	if (scale == "zones") {
		div_input <- div_zones
	} else if (scale == "years_inside") {
		div_input <- div_years[1:28,1:7]
	} else if (scale == "years_outside") {
		div_input <- div_years[30:53,1:7]
	}

	div_out <- iNEXT(div_input, q = c(0, 1, 2), datatype = "abundance")
	asy_est <- as.data.frame(div_out$AsyEst)
	
	# Find the diversity/order column (Order.q or Diversity)
	q_col <- grep("Order.q|Order|order.q|^q$|Diversity", colnames(asy_est), value = TRUE)[1]
	
	evenness_results <- asy_est %>%
		# Filter for Richness and Simpson using both numeric and string labels
		filter(!!sym(q_col) %in% c(0, 2, "Species richness", "Simpson diversity")) %>%
		mutate(q_label = case_when(
			!!sym(q_col) %in% c(0, "Species richness") ~ "q0",
			!!sym(q_col) %in% c(2, "Simpson diversity") ~ "q2"
		)) %>%
		select(Assemblage, q_label, Observed) %>%
		pivot_wider(names_from = q_label, values_from = Observed) %>%
		mutate(Evenness = (q2 - 1) / (q0 - 1)) %>%
		select(Assemblage, q0, q2, Evenness)

	return(evenness_results)
}