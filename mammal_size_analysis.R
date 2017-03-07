# Analysis of individual size distributions of small mammals at NEON

library(dplyr)
library(ggplot2)

load_data = function(path, pattern) {
  files = dir(path, pattern = pattern, full.names = TRUE)
  tables = lapply(files, read.csv, stringsAsFactors = FALSE)
  bind_rows(tables)
}

# NEON.D01.HARV.DP1.10072.001.mam_capturedata.csv manually deleted due to minor data issue
captdata = load_data("data", pattern = "*_capturedata.csv")

ggplot(captdata, aes(x = weight)) +
  geom_histogram() +
  facet_wrap(~siteID, scales = 'free')

# split weight axis into 0.2 ln units and calc total energy use of all indivs in those classes

get_energy_dist <- function(data, max_size) {
  energy_dists = data.frame(siteID = character(0), size_class = numeric(0), energy = numeric(0), stringsAsFactors = FALSE)
  size_class_edges = exp(seq(0, log(max_size) + 0.2, 0.2))
  for (i in seq(1, length(size_class_edges) - 1)){
    class_data = filter(data, weight >= size_class_edges[i] & weight < size_class_edges[i + 1])
    energy = sum(class_data$weight^0.75)
    energy_dists = bind_rows(energy_dists, data.frame(siteID = data$siteID[1],
                                                      size_class = mean(c(log(size_class_edges[i]), log(size_class_edges[i+1]))),
                                                      energy = energy))
  }
  return(energy_dists)
}

max_size = max(captdata$weight, na.rm = TRUE)
energy_dist_data = captdata %>%
  group_by(siteID, nlcdClass) %>%
  do(get_energy_dist(., max_size))

ggplot(energy_dist_data, aes(x = size_class, y = energy)) + 
  geom_bar(stat = 'identity') +
  #facet_grid(siteID~nlcdClass, scales = 'free')
  facet_wrap(~siteID+nlcdClass, scales = 'free')
