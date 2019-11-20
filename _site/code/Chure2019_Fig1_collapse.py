#-*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import mut.viz
import mut.thermo
import matplotlib.pyplot as plt
import seaborn as sns
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
color_palette = sns.color_palette('deep', n_colors=18)
mut.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')

# Note that Razo Mejia et al report the number of tetramers, not dimers
data['repressors'] *= 2 
data = data[data['repressors'] > 0]

# Compute the bohr prameter. 
bohr_wt = mut.thermo.SimpleRepression(R=data['repressors'], 
                                     ep_r=data['binding_energy'],
                                     ka=constants['Ka'], ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'],
                                     effector_conc=data['IPTG_uM']).bohr_parameter()
data['bohr_parameter'] = bohr_wt
grouped = data.groupby(['repressors', 
                        'operator', 
                        'IPTG_uM']).agg(('mean', 'sem')).reset_index()

# Set up the figure canvas. 
fig, ax = plt.subplots(2, 1, figsize=(2, 2.75))


# #############################
# INDUCTION PROFILES 
# #############################
c_range = np.logspace(-2, 4, 500) 
iter = 0
for g, d in grouped.groupby(['repressors', 'operator']):
    ax[0].errorbar(d['IPTG_uM'] / 1E6, d['fold_change_A']['mean'],
                  d['fold_change_A']['sem'], lw=0.25, capsize=0.25,
                  linestyle='none', fmt='.', ms=1, color=color_palette[iter])
    fc = mut.thermo.SimpleRepression(g[0], constants[g[1]],
                                    ka=constants['Ka'], ki=constants['Ki'], 
                                    ep_ai=constants['ep_AI'], 
                                     effector_conc=c_range).fold_change()
    ax[0].plot(c_range / 1E6, fc, lw=0.25, color=color_palette[iter])
    iter += 1
   
    
# #############################
# DATA COLLAPSE
# #############################
iter = 0
for g, d in grouped.groupby(['repressors', 'operator']):
    ax[1].errorbar(d['bohr_parameter']['mean'], d['fold_change_A']['mean'],
                  xerr=d['bohr_parameter']['sem'], yerr=d['fold_change_A']['sem'],
                  lw=0.25, capsize=0.25, linestyle='none', fmt='.', color=color_palette[iter],
                  ms=1)
    iter += 1
bohr_range = np.linspace(-8, 8, 100)
ax[1].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', zorder=1000,
          lw=1)
# #############################
# FORMATTING 
# #############################
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.grid(False)
    a.set_yticks([0, 0.5, 1])
    
ax[0].set_xticks([-6, -4, -2])
ax[1].set_xticks([-5, 0, 5])
ax[0].set_xlabel('IPTG [M]', fontsize=6)
ax[0].set_ylabel('fold-change', fontsize=6)
ax[1].set_xlabel('free energy [$k_BT$]', fontsize=6)
ax[1].set_ylabel('fold-change', fontsize=6)
ax[0].set_xlim([1E-8, 1E-2])
ax[0].set_xscale('log')
ax[1].set_xlim([-8, 8])
plt.tight_layout()
plt.savefig('../../figures/Chure2019_Fig1_collapse.pdf', bbox_inches='tight')


