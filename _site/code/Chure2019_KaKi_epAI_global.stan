/*
* Inference of all allosteric parameters using all induction profiles
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: CC-BY 4.0
*
* Description
* ------------------------------------------------------------------------------
* This model saamples the posterior distribution of the inducer binding
* constants to the active and inactive repressor for a set of J unique inducer
* binding domain mutants. This model assumes that each mutant has a unique
* value for the allosteric energy difference (ep_AI), All other parameter
* values (i.e. DNA binding energy, allosteric energy difference) are taken as
* delta functions at the literature value. 
*/
#include Chure2019_functions.stan
data { 
  // Dimensional parameters
  int<lower=1> N;  
    
  // Architectural parameters
  real<lower=0> R[N]; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  vector<lower=0>[N] c; // Effector concentration.
  real ep_RA[N]; // Binding energy in kBT. 
  int<lower=1> n_sites; // Number of allosteric sites.  

  // Observed parameters.
  vector[N] fc;
  }

parameters {
  real<lower=0> Ka; // log transform of Ka
  real<lower=0> Ki; // log transform of Ki
  real ep_AI;
  real<lower=0> sigma; //  Homoscedastic error
}

transformed parameters {
  real ep_a = log(Ka); 
  real ep_i = log(Ki);
}

model {
  vector[N] mu;
    
  // Define the priors. 
  sigma ~ normal(0, 0.1);
  ep_a ~ normal(2, 2);
  ep_i ~ normal(0, 2);
  ep_AI ~ normal(0, 5);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, ep_RA[i], c[i], ep_a, ep_i, ep_AI, n_sites);
    fc[i] ~ normal(mu[i], sigma);
  }
}