# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import tqdm
constants = mut.thermo.load_constants()

# Load the statistical model. 
model = mut.bayes.StanModel('../stan/Chure2019_empirical_F_inference.stan')

# Define the parameters
n_rep = 10
n_points = 500

# Generate the fake datasets
F_mu = np.linspace(-12, 12, n_points)
fc_mu = (1 + np.exp(-F_mu))**-1

# Draw sigmas out of a half normal 
sig = np.abs(np.random.normal(0, 0.1, len(F_mu)))
dfs = []
for i in range(n_points):
    fc_rand = np.random.normal(fc_mu[i], sig[i], n_rep)
    df = pd.DataFrame([])
    df['fold_change'] = fc_rand
    df['fc_mu'] = fc_mu[i]
    df['fc_sig'] = sig[i]
    df['bohr']  = F_mu[i]
    df['draw'] = i
    dfs.append(df)
df = pd.concat(dfs)

samp_dfs = []
stat_dfs = []
for g, d in tqdm.tqdm(df.groupby(['draw'])):
    # Sample the posterior for the data
    _, samples = model.sample(dict(N=len(d), foldchange=d['fold_change']), 
                            control=dict(adapt_delta=0.95))
    samples['true_mu'] = d['fc_mu'].values[0]
    samples['true_sig'] = d['fc_sig'].values[0]
    samples['true_bohr'] = d['bohr'].values[0]
    samp_dfs.append(samples)

    # Compute the empirical bohr and delta F 
    samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)
    samples['delta_bohr'] = d['bohr'].values[0] - samples['empirical_bohr'] 

    # Compute the delta F error of the reference, given the sigma
    _dbohr_stats = mut.stats.compute_statistics(samples, varnames=['delta_bohr', 
                                                'empirical_bohr', 'fc_mu', 
                                                'fc_sigma'], 
                                                logprob_name='lp__')   
    _dbohr_stats['true_mu'] = d['fc_mu'].values[0] 
    _dbohr_stats['true_sig'] = d['fc_sig'].values[0]
    _dbohr_stats['true_bohr'] = d['bohr'].values[0]
    stat_dfs.append(_dbohr_stats)

samples = pd.concat(samp_dfs)
stats = pd.concat(stat_dfs)
df.to_csv('../../data/Chure2019_empirical_F_simulated_data.csv', index=False)
samples.to_csv('../../data/Chure2019_empirical_F_simulated_data_samples.csv', index=False)
stats.to_csv('../../data/Chure2019_empirical_F_simulated_data_statistics.csv', index=False)