# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.bayes
import mut.thermo
import tqdm
constants = mut.thermo.load_constants()

# Load the data and segregate in to classes
data = pd.read_csv('../../data/Chure2019_compiled_data.csv', comment='#')
Q21M = data[(data['mutant']=='Q21M') & (data['operator']=='O2') &
           (data['repressors']==260)]

# Perform simulations.
n_sim = 800
n_draws = len(Q21M)

# Instantiate a storage data frame
dfs = []
for i in tqdm.tqdm(range(n_sim)):
    # Draw values from the prior
    epRA_draw = np.random.normal(-12,6)
    sigma_draw = np.abs(np.random.normal(0, 0.1))
    
    # Compute the theoretical fold-change. 
    fc_mu = mut.thermo.SimpleRepression(ep_r=epRA_draw, R=260,
                                       ka=constants['Ka'],
                                       ki=constants['Ki'],
                                       effector_conc=Q21M['IPTGuM'], 
                                       n_ns=constants['Nns'],
                                       n_sites=2,
                                       ep_ai=constants['ep_AI']).fold_change()
    
    # Draw N samples from the likelihood 
    fc_draws = np.random.normal(fc_mu, sigma_draw)
    
    # Append to the likelihood
    _df = pd.DataFrame(np.array([fc_draws]).T, columns=['fc_draw'])
    _df['ep_RA'] = epRA_draw
    _df['sigma'] = sigma_draw
    _df['sim_idx'] = i + 1 
    _df['IPTGuM'] = Q21M['IPTGuM'].values
    _df['fc_mu'] = fc_mu.values 
    
    dfs.append(_df)
df = pd.concat(dfs)

# Save the samples to csv. 
df.to_csv('../../data/Chure2019_DNA_prior_predictive_checks.csv',
         index=False)