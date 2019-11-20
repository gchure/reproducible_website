# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import tqdm
import mut.thermo
import mut.bayes
constants = mut.thermo.load_constants()

# Load the raw data
data = pd.read_csv('../../data/Chure2019_compiled_data.csv', comment='#')

# Segregate the data by classifier
DNA_data = data[data['class']=='DNA'].copy()
IND_data = data[data['class']=='IND'].copy()

# Load the Stan models. 
DNA_model = mut.bayes.StanModel('../stan/Chure2019_DNA_binding_energy.stan', 
                                force_compile=True)
KaKi_model = mut.bayes.StanModel('../stan/Chure2019_KaKi_only.stan', 
                                force_compile=True)
KaKi_epAI_model = mut.bayes.StanModel('../stan/Chure2019_KaKi_epAI.stan', 
                                force_compile=True)
empirical_F_model = mut.bayes.StanModel('../stan/Chure2019_empirical_F_inference.stan',
                                force_compile=True)

# ##############################################################################
#  DNA BINDING ENERGY INFERENCE 
# ##############################################################################
mutant_dfs = []
summary_dfs = []
print('Beginning inference of DNA binding energy...')
for g, d in tqdm.tqdm(DNA_data.groupby(['mutant', 'repressors'])):
    # Assemble the data dictionary.
    data_dict = {'J':1,
                 'N': len(d),
                 'idx': np.ones(len(d)).astype(int),
                 'R': d['repressors'], 
                 'Nns': constants['Nns'],
                 'ep_ai': constants['ep_AI'],
                 'Ka': constants['Ka'],
                 'Ki': constants['Ki'], 
                 'n_sites': constants['n_sites'],
                 'c': d['IPTGuM'],
                 'fc': d['fold_change']}
    
    # Sample!
    samples, samples_df = DNA_model.sample(data_dict=data_dict)
    
    # Get the parameter names and rename 
    parnames = samples.unconstrained_param_names()
    new_names = {'{}[{}]'.format(m.split('.')[0], m.split('.')[1]):'{}'.format(m.split('.')[0]) for m in parnames} 
    samples_df.rename(columns=new_names, inplace=True)
    
    # Add identifiers
    samples_df['repressors'] = g[1]
    samples_df['mutant'] = g[0]
    samples_df['operator'] = d['operator'].unique()[0]
    
    # Compute the summarized dataframe
    _df = samples_df[['ep_RA', 'sigma', 'lp__']].copy()
    summary_df = mut.stats.compute_statistics(_df, logprob_name='lp__')
    summary_df['repressors'] = g[1]
    summary_df['mutant'] = g[0]
    summary_df['operator'] = d['operator'].unique()[0]
  
    # Add to storage vector
    mutant_dfs.append(samples_df) 
    summary_dfs.append(summary_df)   
    
# Combine and save to disk    
mutant_df = pd.concat(mutant_dfs, sort=False)
summary_df = pd.concat(summary_dfs, sort=False)
mutant_df.to_csv('../../data/Chure2019_DNA_binding_energy_samples.csv', index=False)
summary_df.to_csv('../../data/Chure2019_DNA_binding_energy_summary.csv', index=False)
print('finished!')

# ##############################################################################
#  KA AND KI INFERENCE
# ##############################################################################
mutant_dfs = []
summary_dfs = []
print('Beginning inference of Ka and Ki...')
for g, d in tqdm.tqdm(IND_data.groupby(['mutant', 'operator'])):
    # Assemble the data dictionary.
    data_dict = {'J':1,
                 'N': len(d),
                 'idx': np.ones(len(d)).astype(int),
                 'R': d['repressors'], 
                 'Nns': constants['Nns'],
                 'ep_AI': constants['ep_AI'],
                 'ep_RA': constants[g[1]],
                 'n_sites': constants['n_sites'],
                 'c': d['IPTGuM'],
                 'fc': d['fold_change']}
    
    # Sample!
    samples, samples_df = KaKi_model.sample(data_dict=data_dict, iter=2000, 
                                            control=dict(adapt_delta=0.99))
    
    # Get the parameter names and rename 
    parnames = samples.unconstrained_param_names()
    new_names = {'{}[{}]'.format(m.split('.')[0], m.split('.')[1]):'{}'.format(m.split('.')[0]) for m in parnames} 
    new_names['Ka[1]'] = 'Ka'
    new_names['Ki[1]'] = 'Ki'
    samples_df.rename(columns=new_names, inplace=True)
    
    # Add identifiers
    samples_df['operator'] = g[1]
    samples_df['repressors'] = d['repressors'].unique()[0]
    samples_df['mutant'] = g[0]
    
    # Compute the summarized dataframe
    _df = samples_df[['Ka', 'Ki', 'sigma', 'lp__']].copy()
    summary_df = mut.stats.compute_statistics(_df, logprob_name='lp__')
    summary_df['repressors'] = d['repressors'].unique()[0]
    summary_df['mutant'] = g[0]
    summary_df['operator'] = g[1]
     
    # Add to storage vector
    mutant_dfs.append(samples_df)  
    summary_dfs.append(summary_df)   
    
# Combine and save to disk    
mutant_df = pd.concat(mutant_dfs, sort=False)
summary_df = pd.concat(summary_dfs, sort=False)
mutant_df.to_csv('../../data/Chure2019_KaKi_only_samples.csv', index=False)
summary_df.to_csv('../../data/Chure2019_KaKi_only_summary.csv', index=False)
print('...finished!')

# ##############################################################################
#  KA, KI, and EpAI INFERENCE
# ##############################################################################
mutant_dfs = []
summary_dfs = []
print('Beginning inference of Ka, Ki, and EpAI...')
for g, d in tqdm.tqdm(IND_data.groupby(['mutant', 'operator'])):
    # Assemble the data dictionary.
    data_dict = {'J':1,
                 'N': len(d),
                 'idx': np.ones(len(d)).astype(int),
                 'R': d['repressors'], 
                 'Nns': constants['Nns'],
                 'ep_RA': constants[g[1]],
                 'n_sites': constants['n_sites'],
                 'c': d['IPTGuM'],
                 'fc': d['fold_change']}
    # Sample!
    samples, samples_df = KaKi_epAI_model.sample(data_dict=data_dict, iter=4000, 
                                       control=dict(adapt_delta=0.995, max_treedepth=11))
    # Get the parameter names and rename 
    parnames = samples.unconstrained_param_names()
    new_names = {'{}[{}]'.format(m.split('.')[0], m.split('.')[1]):'{}'.format(m.split('.')[0]) for m in parnames} 
    new_names['Ka[1]'] = 'Ka'
    new_names['Ki[1]'] = 'Ki'
    samples_df.rename(columns=new_names, inplace=True)
    
    # Add identifiers
    samples_df['operator'] = g[1]
    samples_df['repressors'] = d['repressors'].unique()[0]
    samples_df['mutant'] = g[0]
    _df = samples_df[['Ka', 'Ki', 'ep_AI', 'sigma', 'lp__']].copy()
    summary_df = mut.stats.compute_statistics(_df, logprob_name='lp__')
    summary_df['repressors'] = d['repressors'].unique()[0]

    summary_df['mutant'] = g[0]
    summary_df['operator'] = g[1]
    
    # Add to storage vector
    mutant_dfs.append(samples_df) 
    summary_dfs.append(summary_df)
    
# Combine and save to disk    
mutant_df = pd.concat(mutant_dfs, sort=False)
summary_df = pd.concat(summary_dfs, sort=False)
mutant_df.to_csv('../../data/Chure2019_KaKi_epAI_samples.csv', index=False)
summary_df.to_csv('../../data/Chure2019_KaKi_epAI_summary.csv', index=False)
print('...finished!')

# ##############################################################################
# INFERENCE OF THE EMPIRICAL FREE ENERGY
# ##############################################################################

# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr

# Assign unique identifiers. 
idx = data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM']).ngroup() + 1 
data['idx'] = idx
data.sort_values('idx', inplace=True)
samples_dfs = []

# Make a storage list for the individual statistics
fc_stats = []

# Iter through each grouping and infer
print('Beginning inference of empirical F...')
for g, d in tqdm.tqdm(data.groupby(['mutant', 'repressors', 
                                    'operator', 'IPTGuM'])):
    # Define parameters of the reference state
    ref = d['ref_bohr'].unique()[0]

    # Assemble the data dictionary and sample the posterior
    data_dict = {'N':len(d),
                 'foldchange': d['fold_change']} 
    fit, samples = empirical_F_model.sample(data_dict, iter=5000,
                                            control=dict(adapt_delta=0.99))

    # Compute the empirical bohr and delta F 
    samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)
    samples['delta_bohr'] = samples['empirical_bohr'] - ref
    _dbohr_stats = mut.stats.compute_statistics(samples, 
            varnames=['empirical_bohr', 'fc_mu', 'fc_sigma', 'delta_bohr'],  
            logprob_name='lp__')    
    _dbohr_stats['mutant'] = g[0]
    _dbohr_stats['repressors'] = g[1]
    _dbohr_stats['operator'] = g[2]
    _dbohr_stats['IPTGuM'] = g[3]
    _dbohr_stats['class'] = d['class'].unique()[0]
    fc_stats.append(_dbohr_stats)

# Concatenate and save    
fc_stats = pd.concat(fc_stats)
fc_stats.to_csv('../../data/Chure2019_empirical_F_statistics.csv', index=False)
print('...finished!')

# ##############################################################################
print('All parameter inference has completed. Data files have been saved to ../../data/')
