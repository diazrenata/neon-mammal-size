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
    return(data.frame(size_class = rich_dens$x, rich_prop = rich_dens$y / sum(rich_dens$y)))
  } else {
    return(data.frame(size_class = NULL, prop = NULL))
  }
}

max_size = max(captdata$weight, na.rm = TRUE)
energy_dist_data.2016 = capt.2016 %>%
  group_by(siteID) %>%
  do(get_energy_dist(., max_size)) %>% 
  rename(energy_prop = energy)

rich_dens = captdata %>%
  group_by(siteID) %>%
  do(get_rich_density(.))

energy_rich_data = inner_join(energy_dist_data, rich_dens, by = c("siteID", "size_class"))
  
ggplot(energy_rich_data, aes(x = size_class, y = energy_prop)) +
  geom_bar(stat = 'identity') +
  geom_line(data = rich_dens, mapping = aes(x = size_class, y = rich_prop)) +
  geom_point(data = rich_dens, mapping = aes(x = size_class, y = rich_prop)) +
  #facet_grid(year~siteID, scales = 'free')
  facet_wrap(~siteID, scales = 'free')

ggsave('results/energy_dist_fig.png')

ggplot(energy_rich_data, aes(x = energy_prop, y = rich_prop)) +
  geom_point() +
  facet_wrap(~siteID) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'gam')

ggsave('results/energy_rich_fig.png')

coeff_dets = energy_rich_data %>%
  group_by(siteID) %>%
  summarize(coeff_det = 1 - sum((rich_prop - energy_prop)^2) / sum((rich_prop - mean(rich_prop))^2), corr_coeff = cor(energy_prop, rich_prop))

coeff_det_binned = hist(coeff_dets$coeff_det, breaks = c(min(coeff_dets$coeff_det) - 1, 0, 0.25, 0.5, 0.75, 1))
coeff_det_data = data.frame(breaks = c(-0.125, 0.125, 0.375, 0.625, 0.875), counts = coeff_det_binned$counts)
ggplot(coeff_det_data, aes(x = breaks, y = counts)) +
  geom_bar(stat = 'identity') +
  xlab('Coefficient of Determination') +
  ylab('Number of Species')

ggsave('results/coeff_det_fig.png')


# RMD playing around

# using all years of data, I guess
# what i want is a rank abundance distribution
# abundance measured in biomass and energy

# for each site, for each species, sum all weights

capt.biomass = captdata
capt.biomass = capt.biomass %>% group_by(siteID, taxonID)
capt.biomass = mutate(capt.biomass, sumwt = sum(weight))
capt.biomass = select(capt.biomass, siteID, taxonID, sumwt)
capt.biomass = distinct(capt.biomass)
capt.biomass = capt.biomass %>% ungroup()
capt.biomass = capt.biomass %>% group_by(siteID)
capt.biomass = capt.biomass %>% arrange(desc(sumwt), .by_group = TRUE)
capt.biomass = mutate(capt.biomass, rank = row_number())
capt.biomass = capt.biomass %>% ungroup()
capt.biomass = mutate(capt.biomass, logwt = log(sumwt))

ggplot(capt.biomass, aes(x = rank, y = logwt)) +
  geom_point() +
  facet_wrap(~siteID)


capt.abundance = captdata %>% 
  group_by(siteID, taxonID) %>%
  mutate(abund = n()) %>%
  select(siteID, taxonID, abund)%>%
  distinct() %>%
  ungroup() %>%
  group_by(siteID) %>%
  arrange(desc(abund), .by_group = TRUE) %>%
  mutate(rank.abund = row_number()) %>%
  ungroup() %>% 
  mutate(log.abund = log(abund))



ggplot(capt.abundance, aes(x = rank.abund, y = log.abund)) +
  geom_point() +
  facet_wrap(~siteID)

