# Analysis of individual size distributions of small mammals at NEON

library(dplyr)
library(tidyr)
library(ggplot2)

load_data = function(path, pattern) {
  files = dir(path, pattern = pattern, full.names = TRUE)
  tables = lapply(files, read.csv, stringsAsFactors = FALSE)
  tables_consistent = lapply(tables, select, siteID, date, taxonID, weight)
  bind_rows(tables_consistent)
}

# NEON.D01.HARV.DP1.10072.001.mam_capturedata.csv manually deleted due to minor data issue
# Data from 2015 and 2016 received from NEON. Unzipped and added `_capturedata` to end of
# file name
captdata = load_data("data", pattern = "*_capturedata.csv") %>%
  filter(!is.na(weight))
captdata = separate(captdata, date, c("year", "mo", "day"), sep = "-")

# split weight axis into 0.2 ln units and calc total energy use of all indivs in those classes

get_energy_dist <- function(data, max_size) {
  energy_dists = data.frame(siteID = character(0), size_class = numeric(0), energy = numeric(0), stringsAsFactors = FALSE)
  size_class_edges = exp(seq(0, log(max_size) + 0.2, 0.2))
  total_energy = sum(data$weight^0.75, na.rm = TRUE)
  for (i in seq(1, length(size_class_edges) - 1)){
    class_data = filter(data, weight >= size_class_edges[i] & weight < size_class_edges[i + 1])
    energy = sum(class_data$weight^0.75, na.rm = TRUE) / total_energy
    energy_dists = bind_rows(energy_dists, data.frame(siteID = data$siteID[1],
                                                      size_class = mean(c(log(size_class_edges[i]), log(size_class_edges[i+1]))),
                                                      energy = energy))
  }
  return(energy_dists)
}

get_rich_density <- function(data){
  size_class_edges = seq(0, log(max_size) + 0.2, 0.2)
  size_class_centers = size_class_edges[1:(length(size_class_edges) - 1)] + diff(size_class_edges) / 2
  rich_dens_data = data %>%
    group_by(taxonID) %>%
    filter(weight != 0) %>% # There is least one zero in GUAN
    summarize(mean_mass = mean(weight, na.rm = TRUE),
              mean_log_mass = mean(log(weight), na.rm = TRUE),
              var_mass = var(weight, na.rm = TRUE),
              var_log_mass = var(log(weight), na.rm = TRUE)) %>%
    filter(!is.na(mean_log_mass))
  bw = mean(sqrt(rich_dens_data$var_log_mass), na.rm = TRUE)
  if (nrow(rich_dens_data) > 1){
    rich_dens = density(rich_dens_data$mean_log_mass, bw = bw, n = 28, from = 0.1, to = 5.5)
    return(data.frame(logmass = rich_dens$x, prop = rich_dens$y / sum(rich_dens$y)))
  } else {
    return(data.frame(logmass = NULL, prop = NULL))
  }
}

max_size = max(captdata$weight, na.rm = TRUE)
energy_dist_data = captdata %>%
  group_by(siteID) %>%
  do(get_energy_dist(., max_size)) %>% 
  rename(prop = energy)

rich_dens = captdata %>%
  group_by(siteID) %>%
  do(get_rich_density(.))
  
ggplot(energy_dist_data, aes(x = size_class, y = prop)) + 
  geom_bar(stat = 'identity') +
  geom_line(data = rich_dens, mapping = aes(x = logmass, y = prop)) +
  geom_point(data = rich_dens, mapping = aes(x = logmass, y = prop)) +
  #facet_grid(year~siteID, scales = 'free')
  facet_wrap(~siteID, scales = 'free')
