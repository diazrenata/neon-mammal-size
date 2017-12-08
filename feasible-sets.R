#' RMD 12/8/17
#' Feasible sets for SADs of NEON mammal communities. 

# Load NEON mammal communities...

library(dplyr)
library(tidyr)
library(ggplot2)

# get NEON capture data
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

site_ns = captdata %>%
  group_by(siteID) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  distinct(siteID, n_total)

# pick a site
thisSite= "ABBY"
site = captdata %>% filter(siteID == thisSite)

site_species = site %>% 
  select(taxonID, weight) %>% 
  group_by(taxonID) %>%
  mutate(mean_wt = mean(weight)) %>%
  mutate(n_ind = n()) %>%
  ungroup() %>%
  distinct(taxonID, mean_wt, n_ind) %>%
  mutate(ind_energy = mean_wt ^ 0.75) %>%
  mutate(sp_energy = ind_energy * n_ind) # Note that sp_energy differs from summing all individuals

# site = site %>% mutate(energy = weight ^ 0.75)

# get # of individuals, total energy, and species
site_n = as.numeric(nrow(site))
# site_m = sum(site$energy)
site_m = sum(site_species$sp_energy) # Note that site_m differs from computing directly from capture weight data
site_s = as.numeric(nrow(site_species))

# for site_n, site_s

real_sad = site_species %>%
  arrange(desc(n_ind)) %>%
  pull(n_ind) # don't pull() to keep the df with the rest of the community data...But we will be getting the rest of that later.

library(partitions)
# interlude while I figure out partitions
P(12)
P(5)
parts(5)
diffparts(5)
restrictedparts(5,3,include.zero =F, decreasing = T)

# K, for SADs 
# restrictedparts(N, S, include.zero = F, decreasing = T) returns a matrix of unique RADs for 
# N individuals across S species, with all species represented.

feasible_set = restrictedparts(site_n, site_s, include.zero = F, decreasing = T)

