# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import tqdm
constants = mut.thermo.load_constants()

# Add details to the constants to keep track of new parameters
constants['HG104'] = 22
constants['RBS1L'] = 1740
constants['auto'] = 0
constants['delta'] = 0

# Load the compiled data
flow_data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
flow_data['repressors'] *= 2
flow_data.rename(columns={'IPTG_uM':'IPTGuM', 'fold_change_A':'fold_change'},
                inplace=True)
flow_data['method'] = 'flow cytometry'
flow_data = flow_data[flow_data['repressors'] > 0].copy()

 ##############################################################################
# INFERENCE
# ##############################################################################
# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr

# Assign unique identifiers. 
idx = data.groupby(['method', 'repressors', 'operator', 'IPTGuM']).ngroup() + 1 
data['idx'] = idx
data.sort_values('idx', inplace=True)
samples_dfs = []
def infer_empirical_bohr(data, model, groupby=['repressors', 'operator', 'IPTGuM', 'method'],
                        verbose=True, force_compile=False, **kwargs):
    """
    Infers the empirical bohr parameter (and relevant correction) for a collection of 
    fold-change measurements
    
    Parameters
    ----------
    data: pandas DataFrame object
        The data from which the empirical bohr will be determined. This should have at least
        a fold-change column and a grouping parameter.
    model: str
        Path to Stan model to load. Model will be compiled if `force_compile`==True.
    groupby: list, optional
        List of identifiers by which to group the supplied data. Default groups by 
        'mutant', 'repressors', 'operator', and 'IPTGuM'
    verbose: bool
        If true, the progress will be printed to screen as a bar. 
    force_compile: bool
        If True, the stan model will be recompiled.
    **kwargs: keyword arguments
        kwargs to be passed to the sampler.
        
    Returns
    -------
    statistics: pandas DataFrame
        Dataframe of statistics for relevant parameters.
    """
    
    # Load the stan model and compile if needed. 
    model = mut.bayes.StanModel(model, force_compile=force_compile)
    
    # Make a storage list for the individual statistics
    fc_stats = []
    
    # Make a quiet or loud iterator. 
    if verbose:
        iter = tqdm.tqdm(data.groupby(groupby))
    else:
        iter = data.groupby(groupby)
        
    # Iter through each grouping and infer
    for g, d in iter: 
        # Define parameters of the reference state
        ref = d['ref_bohr'].unique()[0]
    
        # Assemble the data dictionary and sample the posterior
        data_dict = {'N':len(d),
                     'foldchange': d['fold_change']} 
        fit, samples = model.sample(data_dict, **kwargs)

        # Compute the empirical bohr and delta F 
        samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)

        # Identify the extrema
        samples['delta_bohr'] = samples['empirical_bohr'] - ref

        _dbohr_stats = mut.stats.compute_statistics(samples, 
                varnames=['empirical_bohr', 'fc_mu', 'fc_sigma', 'delta_bohr'],  
                logprob_name='lp__')    
        _dbohr_stats['method'] = g[-1]
        _dbohr_stats['repressors'] = g[0]
        _dbohr_stats['operator'] = g[1]
        _dbohr_stats['IPTGuM'] = g[2]
        _dbohr_stats['pred_bohr'] = ref
        fc_stats.append(_dbohr_stats)
    
    return pd.concat(fc_stats)


fc_stats = infer_empirical_bohr(data[data['method']=='flow cytometry'], '../stan/Chure2019_empirical_F_inference.stan', 
            **dict(iter=5000, control=dict(adapt_delta=0.99)))
fc_stats.to_csv('../../data/csv/RazoMejia2018_empirical_F_statistics.csv', index=False)


