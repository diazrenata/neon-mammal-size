# RMD December 2017
# simulating potential communities using the sizes of the species found at the site

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

# pick a site
thisSite= "CPER"
site = captdata %>% filter(siteID == thisSite)

site_species = site %>% 
  select(taxonID, weight) %>% 
  group_by(taxonID) %>%
  mutate(mean_wt = mean(weight)) %>%
  ungroup() %>%
  distinct(taxonID, mean_wt) %>%
  mutate(ind_energy = mean_wt ^ 0.75)

site = site %>% mutate(energy = weight ^ 0.75)

# get # of individuals, total energy, and species
site_n = as.numeric(nrow(site))
site_m = sum(site$energy)
site_s = as.numeric(nrow(site_species))


# # I think the species-body size distribution will be important....
# hist(site_species$mean_wt )



# try and generate just one fake community, then look into iterations & saving successful runs.


kept_comm = list()
comm_specs = matrix(nrow = 0, ncol = 3)
colnames(comm_specs) =  c("n", "m", "s")
set.seed(10)

#while (keep == FALSE) {

i = 1

while(i < 1000) {  

    sim_comm = site_species %>%
    mutate(nind = 0) %>%
    mutate(sp_energy = 0)
   
    sim_n = 0
    sim_m = 0
    sim_s = 0
    
   # while(sim_n <= site_n && sim_m <= site_m) {
    while(sim_m < site_m) {
    # new individuals belong to any of the 15 species
    new.ind = as.numeric(sample(site_s, size = 1, replace = TRUE, prob = NULL))
    sim_comm[new.ind, 4] <- sim_comm[new.ind, 4] + 1
    sim_comm[new.ind, 5] <- sim_comm[new.ind, 4] * sim_comm[new.ind, 3]
    sim_n = sum(sim_comm[,4])
    sim_m = sum(sim_comm[,5])
    sim_s = length(which(sim_comm[,4] > 0))
  }
  comm_specs = rbind(comm_specs, c(sim_n, sim_m, sim_s))
  
  #if (abs(site_n - sim_n) <= (.2*site_n) && abs(site_m - sim_m) <= .2*site_m && min(sim_comm[,4]) >0) {
    kept_comm[[length(kept_comm) + 1]] <- sim_comm
 # } 
  i = i + 1
}
comm_specs <- as.data.frame(comm_specs)
comm_specs <- cbind(comm_specs, abs(site_n - comm_specs[,1]))

comm_specs <- cbind(comm_specs, abs(site_m - comm_specs[,2]))
comm_specs <- cbind(comm_specs, abs(site_s - comm_specs[,3]))
colnames(comm_specs)[4:6] <- c('diff_n', 'diff_m', 'diff_s')

detach("package:plyr", unload=TRUE) 

site_abund = site %>% 
  group_by(taxonID) %>%
  mutate(nind = n()) %>%
  select(taxonID, nind)%>%
  distinct() %>%
  ungroup() %>%
  arrange(desc(nind)) %>%
  mutate(rank_abund = row_number())

site_energy = site %>% 
  group_by(taxonID) %>%
  mutate(sp_energy = sum(energy)) %>%
  select(taxonID, sp_energy)%>%
  distinct() %>%
  ungroup() %>%
  arrange(desc(sp_energy)) %>%
  mutate(rank_energy = row_number())


# library(plyr)
ranks = cbind(site_energy[,2:3], site_abund[,2:3])
              
ranks = as_tibble(ranks) %>%
  select(rank_energy, sp_energy, rank_abund, nind) %>%
  mutate(simulation = 'real')



keep = comm_specs %>%
  mutate(iteration = row_number()) %>%
  filter(diff_n == 1, diff_s == 0)
# keep = comm_specs %>%
#   mutate(iteration = row_number()) %>%
#   filter(diff_n < 470, diff_m <5, diff_s == 0)

plot_comms = kept_comm[keep$iteration]


for(i in 1:10) {
  thisone = plot_comms[[i]] 
  thisone = thisone %>% 
    arrange(desc(sp_energy)) %>%
    mutate(rank_energy = 1:nrow(thisone)) %>%
    arrange(desc(nind)) %>%
    mutate(rank_abund = 1:nrow(thisone))
  plot_comms[[i]] <- thisone
  
  thisone = thisone %>%
    select(rank_energy, sp_energy, rank_abund, nind) %>%
    mutate(simulation = paste0('sim_', i))
  
  ranks = rbind(ranks, thisone)
}

ranks = ranks %>% 
  mutate(log_ind = log(nind)) %>%
  mutate(log_energy = log(sp_energy))

ggplot(ranks, aes(x = rank_abund, y = log_ind)) +
  geom_point() +
  facet_wrap(~simulation)


ggplot(ranks, aes(x = rank_energy, y = log_energy)) +
  geom_point() +
  facet_wrap(~simulation)

