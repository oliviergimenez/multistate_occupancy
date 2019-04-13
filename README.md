# Fitting a multistate occupancy model with uncertainty in Jags

## Motivation

We would like to implement a multistate single-season occupancy model (Nichols et al. 2007) in Jags. To do so, we adopt a hidden Markov modeling formulation of the model (Gimenez et al. 2014 and associated [Wiki](http://occupancyinesurge.wikidot.com/multiple-states-and-uncertainty)). To illustrate the analysis, we use simulated data provided by Donovan et al. (2007).

## Data

The states are 1 for site unoccupied, 2 for occupied with no production of young ('uncertain' non-breeding state) and 3 for occupied with successful reproduction ('certain' breeding state). The observations are 0 for species not observed, 1 for species observed and 2 for species observed with young.

These data were simulated with the following parameter values:
* probability that the site is occupied by non-breeders psi1 = 0.3;
* probability that the site is occupied by breeders psi2 = 0.5;
* detection probability of non-breeders p1 = 0.5;
* detection probability of breeders p2 = 0.7;
* probability of detecting evidence of reproduction, given the site is occupied with young delta = 0.8;
* number of sites R = 250.

```{r message=FALSE, warning=FALSE, paged.print=FALSE}
dat <- readr::read_csv2('https://raw.githubusercontent.com/oliviergimenez/multistate_occupancy/master/multiocc.csv')
#head(dat)
#tail(dat)
#sum(dat==1)
#sum(dat==2)
#sum(dat==0)
```
