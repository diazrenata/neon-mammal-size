# NEON small mammal community simulations
RMD 12.7.17
## General idea
* Trying to look at quasi feasible sets of rank-energy and rank-abundance distributions, working with NEON small mammal data.
* Do simulated communities with the same N, S, and E as real communities have comparable rank-abundance and rank-energy distributions?
* Inspired by a feasible set approach, and trying to deal with the conflict between combinatorics and continuous data (biomass, energy use)

## Approach
* For potentially many NEON sites (code is general)...
  * Pull N (no. of individuals), S (no. of species), E (total community metabolic rate) for the whole community, and mean body size and metabolic rate for each species at that sites
    * energy use = weight ^ .75
  * Generate many simulated randomly assembled communities and select the ones with similar S, N, E to real
    * Each individual assigned randomly to a species until total N or E = real N or E
    * Select the communities with N = real N +/- 5, E = real E +/-5, S = real S
    * Plot rank-abundance and rank-energy distributions for selected simulated communities and compare to real RAD, REDs.
  * Out of curiosity I also started removing the E constraint.

## Results so far
* 
