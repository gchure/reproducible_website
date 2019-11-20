# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the prior predictive checks.
data = pd.read_csv('../../data/Chure2019_empirical_F_prior_predictive_checks.csv')

# Define the true distributions
x = np.linspace(0, 1, 100)
fc_mu = np.ones(len(x)) * 1 / len(x)
fc_sigma = np.abs(scipy.stats.norm(0, 0.1).pdf(x))

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
ax = [ax1, ax2, ax3]
for a in ax:
   a.xaxis.set_tick_params(labelsize=6)
   a.yaxis.set_tick_params(labelsize=6)
ax3.set_xlim(-0.3, 1.3)
# Adjust limits
ax1.set_ylim([0, 0.0155])
ax2.set_ylim([0, 6])
ax3.set_ylim([0, 1])
ax1.set_xlim([0, 1])
ax2.set_xlim([0, 0.5])
ax1.set_yticks([])
ax2.set_yticks([])

# Add labels
ax1.set_xlabel('$\mu$', fontsize=8)
ax2.set_xlabel('$\sigma$', fontsize=8)
ax3.set_xlabel('fold-change', fontsize=8)
ax3.set_ylabel('cumulative distribution', fontsize=8)

# Add panel labels. 
fig.text(-0.05, 0.95, '(A)', fontsize=8)
fig.text(-0.05, 0.5, '(B)', fontsize=8)
# ##############################################################################
# PRIOR DISTRIBUTIONS
# ##############################################################################
ax1.plot(x, fc_mu, '-', color=colors['red'])
ax2.plot(x, fc_sigma, '-', color=colors['red'])
ax2.fill_between(x, fc_sigma, '-', color=colors['light_red'])
ax1.fill_between(x, fc_mu, color=colors['light_red'])

# ##############################################################################
# PRIOR SAMPLES 
# ##############################################################################
x1 = np.random.normal(0, 0.0005, len(data['mu'].unique())) + 0.014
x2 = np.random.normal(0, 0.17, len(data['mu'].unique())) + 5.5
ax1.hlines(0.014, 0, 1, lw=35, color='w')
ax1.plot(data['mu'].unique(), x1, 'k.', ms=0.8)
ax2.plot(data['sigma'].unique(), x2, 'k.', ms=0.8)
ax2.hlines(5.5, 0, 1, lw=35, color='w')

# ##############################################################################
# PERCENTILES
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
perc_cmap = sns.color_palette('Reds', n_colors=(len(percs) + 2))
zorder = [11, 12, 13, 14, 15, 16, 17]
z = {p:z for p, z in zip(percs, zorder)}
c = {p:c for p, c in zip(percs, perc_cmap)}

# Group by the ECDFS
dfs = []
for g, d in data.groupby(['draw']):
    d = d.copy()
    d = d.sort_values('fold_change')
    d['y'] = np.arange(0, 10, 1) / 10
    dfs.append(d)
_data = pd.concat(dfs)


df = pd.DataFrame([])
for g, d in _data.groupby('y'):
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(data['fold_change'], [low, upper])
        df = df.append({'percentile': p,
                    'y': g + 0.1,
                   'fc_low':_percs[0],
                   'fc_high': _percs[1]},
                   ignore_index=True)

for g, d in df.groupby('percentile'):
   ax3.fill_betweenx(np.linspace(0, 1.1, 10), d['fc_low'], d['fc_high'], 
                color=c[g], zorder=z[g], label=int(g))
leg = ax3.legend(title='percentile', fontsize=8)
leg.set_zorder(1001)
leg.get_title().set_fontsize(8)
plt.tight_layout()
plt.savefig('../../figures/Chure2019_FigS7_empirical_F_prior_predictive_checks.pdf',
             bbox_inches='tight')


