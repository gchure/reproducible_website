# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import mut.stats
constants = mut.thermo.load_constants()

# Load the data and define the stan models.
data = pd.read_csv('../../data/Chure2019_compiled_data.csv', comment='#')
epRA_model = mut.bayes.StanModel('../stan/Chure2019_DNA_binding_energy_global.stan')
allo_model = mut.bayes.StanModel('../stan/Chure2019_KaKi_epAI_global.stan')


# ##############################################################################
# INFERENCE OF DNA BINDING ENERGY
# ##############################################################################
DNA_data = data[data['class']=='DNA']
samples_dfs = []
stats_dfs = []
for g, d in DNA_data.groupby(['mutant']):
    d = d.copy()

    # Define the data dictionary. 
    data_dict = {'N':len(d), 
                 'R':d['repressors'],
                 'ep_ai':constants['ep_AI'],
                 'Ka':constants['Ka'],
                 'Ki':constants['Ki'],
                 'Nns':constants['Nns'],
                 'n_sites': constants['n_sites'],
                 'c':d['IPTGuM'],
                 'fc':d['fold_change']}

    # Sample!
    fit, samples = epRA_model.sample(data_dict)
    stats = mut.stats.compute_statistics(samples, varnames=['sigma', 'ep_RA'],
                                         logprob_name='lp__')
    samples['mutant'] = g
    stats['mutant'] = g
    samples_dfs.append(samples)
    stats_dfs.append(stats)

# Concatenate the data frames and save.
samples = pd.concat(samples_dfs)
stats = pd.concat(stats_dfs)
samples.to_csv('../../data/Chure2019_global_DNA_binding_energy_samples.csv',
               index=False)
stats.to_csv('../../data/Chure2019_global_DNA_binding_energy_summary.csv',
               index=False)

# ##############################################################################
# ALLOSTERIC PARAMETER SAMPLING
# ##############################################################################
IND_data = data[data['class']=='IND']
samples_dfs = []
stats_dfs = []
for g, d in IND_data.groupby(['mutant']):
    d = d.copy()
    ops = [constants[o] for o in d['operator'].values]
    # Define the data dictionary. 
    data_dict = {'N':len(d), 
                 'R':d['repressors'],
                 'ep_RA': ops,
                 'Nns':constants['Nns'],
                 'n_sites': constants['n_sites'],
                 'c':d['IPTGuM'],
                 'fc':d['fold_change']}

    # Sample!
    fit, samples = allo_model.sample(data_dict, control=dict(adapt_delta=0.95))
    stats = mut.stats.compute_statistics(samples, varnames=['Ka', 'Ki', 
                                                  'ep_AI', 'sigma'],
                                                  logprob_name='lp__')
    samples['mutant'] = g
    stats['mutant'] = g
    samples_dfs.append(samples)
    stats_dfs.append(stats)

# Concatenate the data frames and save.
samples = pd.concat(samples_dfs)
stats = pd.concat(stats_dfs)
samples.to_csv('../../data/Chure2019_global_KaKi_epAI_samples.csv',
               index=False)
stats.to_csv('../../data/Chure2019_global_KaKi_epAI_summary.csv',
               index=False)

