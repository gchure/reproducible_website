# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import mut.stats
import joblib
import multiprocessing as mp
cpus = mp.cpu_count() - 2
import tqdm
constants = mut.thermo.load_constants()

# Load the prior predictive check data. 
prior_data = pd.read_csv('../../data/Chure2019_IND_prior_predictive_checks.csv')

# Load the stan model. 
KaKi_model = mut.bayes.StanModel('../stan/Chure2019_KaKi_only.stan')
KaKi_epAI_model = mut.bayes.StanModel('../stan/Chure2019_KaKi_epAI.stan')
_model = {'KaKi_only': KaKi_model, 'KaKi_epAI':KaKi_epAI_model}

# Set up a dataframe to store the properties.
samples_dfs = []
sbc_dfs = []

# Define the thinning constant for computing the rank statistic. 
thin = 5

# Set up a single round of the SBC as a function for easy parallelization
def sbc(g, d): 
        # Generate the data dictionary. 
        data_dict = {'J':1,
                    'N': len(d),
                    'idx': np.ones(len(d)).astype(int),
                    'ep_RA': -13.9,
                    'R': np.ones(len(d)) * constants['RBS1027'],
                    'Nns': 4.6E6,
                    'n_sites': constants['n_sites'],
                    'c': d['IPTGuM'],
                    'fc': d['fc_draw']}
        # Define the columns for renaming
        columns={'Ka[1]': 'Ka', 'sigma[1]':'sigma', 'Ki[1]':'Ki',
                 'ep_a[1]':'ep_a', 'ep_i[1]': 'ep_i'}

    
        # Determine the ground truth for each parameter.
        gt = {'Ka': d['ka'].unique(),
                  'Ki': d['ki'].unique(),
                  'ep_a':d['ep_a'].unique(),
                  'ep_i':d['ep_i'].unique(),
                  'sigma': d['sigma'].unique()}
        if g[0] == 'KaKi_only': 
            data_dict['ep_AI'] = constants['ep_AI']
            pars = ['Ka', 'Ki', 'ep_a', 'ep_i', 'sigma']
        else:
            gt['ep_AI'] = d['ep_ai'].unique()
    
            pars = ['Ka', 'Ki', 'ep_AI', 'ep_a', 'ep_i', 'sigma']
            columns['ep_AI[1]'] = 'ep_AI'
    
            
        # Sample the model
        model = _model[g[0]]
        _, samples = model.sample(data_dict=data_dict,  iter=2000, n_jobs=1, chains=4)
        samples.rename(columns=columns, inplace=True)
        samples['sim_idx'] = g[1]
        samples['model'] = g[0]
        samples_dfs.append(samples)
        
        # Compute the properties for each parameter. 
        _sbc_dfs = []
        for p in pars:
            print(p, samples[p].head())
            _df = pd.DataFrame([])
            z_score = (np.mean(samples[p]) - gt[p]) / np.std(samples[p])
            shrinkage = 1 - (np.var(samples[p]) / np.var(prior_data[p.lower()].unique()))
            _df['z_score'] = z_score
            _df['shrinkage'] = shrinkage
            _df['param'] = p 
            _df['rank'] = np.sum(samples[p].values[::thin] < gt[p])
            _df['rank_ndraws'] = len(samples[p].values[::thin])
            _df['post_median'] = np.mean(samples[p])
            _df['post_mean'] = np.median(samples[p])
            _df['post_mode'] = samples.iloc[np.argmax(samples['lp__'].values)][p]
            _df['model'] = g[0]
            _df['ground_truth'] = gt[p]
            _sbc_dfs.append(_df)
            
        _sbc_dfs = pd.concat(_sbc_dfs)
        _sbc_dfs['sim_idx'] = g[1]
        sbc_dfs.append(_sbc_dfs) 
        return [samples, _sbc_dfs]

out = joblib.Parallel(n_jobs=cpus)(joblib.delayed(sbc)(g, d)for g, d in tqdm.tqdm(prior_data.groupby(['model', 'draw'])))
_samples = [a[0] for a in out]
_sbc = [a[1] for a in out]
sbc_df = pd.concat(_sbc) 
sbc_df.to_csv('../../data/Chure2019_IND_sbc_samples.csv', index=False)

