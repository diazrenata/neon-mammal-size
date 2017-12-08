# NEON small mammal community simulations
RMD 12.7.17
## General idea
* Trying to look at quasi feasible sets of rank-energy and rank-abundance distributions, working with NEON small mammal data.
* Do simulated communities with the same N, S, and E as real communities have comparable rank-abundance and rank-energy distributions?
* Inspired by a feasible set approach, and trying to deal with the conflict between combinatorics and continuous data (biomass, energy use)
  * **this isn't really feasible sets!**
  


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
* Caveats...so far few iterations, only a couple of species-poor communities, having to allow large deviations from real just to debug the code, computationally slow
* Simulated RADs and REDs much flatter than real.
  * RAD seems especially weird to me - prev. feasible sets (RAD given S, N) would indicate otherwise
  * Maybe it's because of the E constraint. Removing it to see what happens
  * Or bc this is  ~9 runs/small communities.

## Other thoughts
* Need to consider intraspecific variation in body size - right now I'm  taking the mean as representative. Should really be pulling from a distribution, but not sure how much I expect this to impact outcomes
* Curious about doing this from a metacommunity/local community perspective.
  * I.e. random pulls from a broader species pool; how many ways can you even maintain E
  * Related to DS/PB/PP histories...
* And datasets for larger communities....although if less species rich communities (like mammals vs. plankton or microbes or even trees) behave qualitatively differently than more species rich ones, *that's* interesting
