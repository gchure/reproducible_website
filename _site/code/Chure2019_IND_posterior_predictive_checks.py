# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import tqdm
import mut.thermo
constants = mut.thermo.load_constants()

# Load the data
data = pd.read_csv('../../data/Chure2019_compiled_data.csv', comment='#')
data = data[(data['mutant']=='Q294K') & (data['operator']=='O2')].copy()
IPTGuM = data['IPTGuM'].unique()

# Load the samples
kaki_only = pd.read_csv('../../data/Chure2019_KaKi_only_samples.csv')
kaki_only = kaki_only[(kaki_only['mutant']=='Q294K') & (kaki_only['operator']=='O2')]
kaki_only['ep_AI'] = constants['ep_AI']
kaki_epAI = pd.read_csv('../../data/Chure2019_KaKi_epAI_samples.csv')
kaki_epAI = kaki_epAI[(kaki_epAI['mutant']=='Q294K') & (kaki_epAI['operator']=='O2')]

# Define the model samples
samples = {'KaKi_only':kaki_only, 'KaKi_epAI':kaki_epAI}

dfs = [] 
for m, s in samples.items():
    print(f'Processing model: {m}')
    counts = data.groupby('IPTGuM').count()['fold_change'].values.astype(int)

    sigma = s['sigma'].values
    for i in tqdm.tqdm(range(len(s))):
        for j in range(len(IPTGuM)):
            fc_mu = mut.thermo.SimpleRepression(R=260, ep_r=-13.9,
                                    ka=s['Ka'].values[i], ki=s['Ki'].values[i],
                                    effector_conc=IPTGuM[j], ep_ai=s['ep_AI'].values[i]).fold_change()

            draws = np.random.normal(fc_mu, sigma[i], size=counts[j])
            _df = pd.DataFrame([], columns = ['IPTGuM'])
            _df['fc_mu'] = draws
            _df['IPTGuM'] = IPTGuM[j]
            _df['samp'] = i
            _df['model'] = m
            dfs.append(_df)    
ppc_df = pd.concat(dfs)
ppc_df.to_csv('../../data/csv/Chure2019_IND_posterior_predictive_samples.csv')
