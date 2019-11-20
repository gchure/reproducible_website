# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.stats
import mut.viz
color = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

#  Load the data
data = pd.read_csv('../../data/Chure2019_summarized_data.csv', comment='#')
data = data[data['class'] == 'IND']
KaKi_only_samples = pd.read_csv('../../data/Chure2019_KaKi_only_samples.csv')
KaKi_epAI_samples = pd.read_csv('../../data/Chure2019_KaKi_epAI_samples.csv')

# Determine the unique repressor copy numbers
ops = np.sort(data['operator'].unique())
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0

# Change this parameter to decide which plot to make
MODEL = 'KaKi_epAI'

# ##############################################################################
#  FIGURE WITH KAKI FIT ONLY
# ##############################################################################
fig, ax = plt.subplots(len(ops),len(ops), figsize=(7,5), sharex=True, sharey=True) 

# Format the axes
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8) 
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xscale('symlog', linthreshx=1E-2)


# Add appropriate labels
for i in range(len(ops)):
    ax[i, 0].set_ylabel('fold-change', fontsize=8)
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=8)
    ax[i, 0].text(-0.5, 0.55, ops[i], fontsize=8, backgroundcolor=pboc['pale_yellow'],
                 transform=ax[i,0].transAxes, rotation='vertical')
    ax[0, i].set_title(ops[i], fontsize=8, backgroundcolor=pboc['pale_yellow'], y=1.08)
    ax[i, i].set_facecolor('#e4e7ec')
for i in range(3):
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    
# Add predictor titles
fig.text(-0.04, 0.53, 'fit strain', fontsize=8, backgroundcolor='#E3DCD0', rotation='vertical')
fig.text(0.435, 0.98, 'comparison strain', fontsize=8, backgroundcolor='#E3DCD0', rotation='horizontal')
    
# Plot the data. 
for g, d in data.groupby(['mutant']):
    g = g.upper()
    for i, _ in enumerate(ops):
        for j, _ in enumerate(ops):
            _d = d[d['operator'] == ops[j]]
            if i == j:
                face = 'w'
                edge = color[g]
            else:
                face = color[g]
                edge = color[g]
            _ = ax[i, j].errorbar(_d['IPTGuM'], _d['mean'], _d['sem'], markerfacecolor=face,
                                 markeredgecolor=edge, color=edge, lw=0.15, linestyle='none', fmt='o',
                                 ms=2.5, label=g)
           
            # Plot the best-fit lines. 
            for k, m in enumerate(data['mutant'].unique()):
                _d = data[(data['mutant']==m) & (data['operator'] == ops[i])]
                
                # Get the binding energies. 
                if MODEL == 'KaKi_only':
                    _samps = KaKi_only_samples[(KaKi_only_samples['mutant']==m) &\
                         (KaKi_only_samples['operator']==ops[i])]
                    _samps['ep_AI'] = 4.5
                else:
                    _samps = KaKi_epAI_samples[(KaKi_epAI_samples['mutant']==m) &\
                        (KaKi_epAI_samples['operator']==ops[i])]
                Ka = _samps['Ka']
                Ki = _samps['Ki']
                epAI = _samps['ep_AI']
                cred_region = np.zeros((2, len(c_range)))
                for z, c in enumerate(c_range): 
                    # Compute the fold-change   
                    fc = mut.thermo.SimpleRepression(R=constants['RBS1027'], 
                                            ep_r=constants[ops[j]], ka=Ka, 
                                                 ki=Ki, ep_ai=epAI,
                                                effector_conc=c, n_sites=constants['n_sites'],
                                                n_ns=constants['Nns']).fold_change()
                    cred_region[:, z] = mut.stats.compute_hpd(fc, 0.95)
                
                # Plot the fit. 
                # _ = ax[i, j].plot(c_range / 1E6, fc[:, 0], color=color[m], lw=0.75) 
                _ = ax[i, j].fill_between(c_range, cred_region[0, :], 
                                   cred_region[1, :], color=color[m], alpha=0.2) 

cred_region[0, :]
cred_region[1, :]
_  = ax[0, 2].legend(fontsize=8, bbox_to_anchor=(1.04, 0.95))
plt.subplots_adjust(wspace=0.05, hspace=0.05)

if MODEL == 'KaKi_only':
   plt.savefig('../../figures/Chure2019_FigS17_KaKi_IND_pairwise_predictions.pdf', 
               bbox_inches='tight')
elif MODEL == 'KaKi_epAI':
   plt.savefig('../../figures/Chure2019_FigS18_KaKi_epAI_IND_pairwise_predictions.pdf', 
            bbox_inches='tight')
    