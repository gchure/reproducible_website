# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.bayes

# Load the data
data = pd.read_csv('../../data/Chure2019_empirical_F_prior_predictive_checks.csv')
data.rename(columns={'mu':'fc_mu', 'sigma':'fc_sigma'}, inplace=True)
model = mut.bayes.StanModel('../stan/Chure2019_empirical_F_inference.stan', 
                            force_compile=True)

# Iterate through each simulation
thin = 5
sbc_dfs = []
samples_dfs = []
for g, d in tqdm.tqdm(data.groupby('draw')):
    
    # Determine the ground truth for each parameter.
    gt = {'fc_mu': d['fc_mu'].unique(),
         'fc_sigma': d['fc_sigma'].unique()}

    # Generate the data dictionary. 
    data_dict = {'N': len(d),
                'foldchange':d['fold_change']}

    # Sample the model
    _, samples = model.sample(data_dict=data_dict)
    samples['sim_idx'] = g
    samples_dfs.append(samples)
    
    # Compute the properties for each parameter. 
    _sbc_dfs = []
    for p in ['fc_mu', 'fc_sigma']:
        _df = pd.DataFrame([])
        z_score = (np.mean(samples[p]) - gt[p]) / np.std(samples[p])
        shrinkage = 1 - (np.var(samples[p]) / np.var(data[p].unique()))
        _df['z_score'] = z_score
        _df['shrinkage'] = shrinkage
        _df['param'] = p 
        _df['rank'] = np.sum(samples[p].values[::thin] < gt[p])
        _df['rank_ndraws'] = len(samples[p].values[::thin])
        _df['post_mean'] = np.mean(samples[p])
        _df['post_median'] = np.median(samples[p])
        _df['post_mode'] = samples.iloc[np.argmax(samples['lp__'].values)][p]
        _df['ground_truth'] = gt[p]
        _sbc_dfs.append(_df)
        
    _sbc_dfs = pd.concat(_sbc_dfs)
    _sbc_dfs['sim_idx'] = g
    sbc_dfs.append(_sbc_dfs) 
sbc_df = pd.concat(sbc_dfs) 
    
sbc_df.to_csv('../../data/Chure2019_empirical_F_sbc_statistics.csv', index=False)


