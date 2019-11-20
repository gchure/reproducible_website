# -*- coding: utf-8 -*- 
import pandas as pd
import numpy as np
import mut.thermo
import mut.stats
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load data to get IPTG concentrations
data = pd.read_csv('../../data/Chure2019_compiled_data.csv')
IPTG = data[(data['mutant']=='Q294K') & (data['operator']=='O2')]['IPTGuM']

# Define the constants relative for drawing samples
n_draws = 1000
ka = np.random.lognormal(2, sigma=2, size=n_draws)
ki = np.random.lognormal(0, sigma=2, size=n_draws)
ep_ai = np.random.normal(0, 5, n_draws)

sigma = np.abs(np.random.normal(0, 0.1, n_draws))
model_names = ['KaKi_only', 'KaKi_epAI']
dfs = []

for m in model_names:
   for i in range(n_draws): 

        if m == 'KaKi_only': 
           args = dict(ka=ka[i], ki=ki[i], ep_ai=constants['ep_AI'])
        else: 
            args = dict(ka=ka[i], ki=ki[i], ep_ai=ep_ai[i])
        arch = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, effector_conc=IPTG,
                                           **args).fold_change() 
        _df = pd.DataFrame([]) 
        _df['fc_draw'] = np.random.normal(arch, sigma[i])
        _df['ep_a'] = np.log(ka[i])
        _df['ka'] = ka[i]
        _df['ep_i'] = np.log(ki[i])
        _df['ki'] = ki[i]
        _df['ep_ai'] = args['ep_ai'] 
        _df['sigma'] = sigma[i]
        _df['draw'] = i
        _df['model'] = m
        _df['IPTGuM'] = IPTG.values
        dfs.append(_df)
df = pd.concat(dfs)
df.to_csv('../../data/Chure2019_IND_prior_predictive_checks.csv', index=False)