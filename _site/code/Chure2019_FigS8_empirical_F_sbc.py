# -*- coding: uttf-8 -*-
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the SBC data
data = pd.read_csv('../../data/Chure2019_empirical_F_sbc_samples.csv')

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# Set limits
ax[0].set_xlim([0, 1.05])
ax[0].set_ylim([-5, 5])
ax[1].set_xlim([0, 800])
ax[1].set_ylim([0, 1])

# Set labels
ax[0].set_xlabel('shrinkage', fontsize=8)
ax[0].set_ylabel('z-score', fontsize=8)
ax[1].set_xlabel('rank statistic', fontsize=8)
ax[1].set_ylabel('cumulative distribution', fontsize=8)

# Add panel labels
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# ##############################################################################
# SENSITIVITY
# ##############################################################################
param_colors = {'fc_mu':colors['red'], 'fc_sigma': 'k'}
legend = {'fc_mu': 'µ', 'fc_sigma': '$\sigma$'}
for g, d in data.groupby('param'):
    if g == 'fc_sigma':
        z = 100
    else:
        z = 101
    ax[0].plot(d['shrinkage'], d['z_score'], '.', color=param_colors[g], 
                ms=1, zorder=z, label=legend[g])

leg = ax[0].legend(title='parameter', fontsize=8)
leg.get_title().set_fontsize(8)

# ##############################################################################
# TRUE UNIFORM DISTRIBUTION
# ##############################################################################
n_sim = data.sim_idx.max()
L = np.arange(0, n_sim, 1)
R = data['rank_ndraws'].unique() 

# Envelope of cdf 99%
y = scipy.stats.randint.cdf(L, 0, R)
std = np.sqrt(y * (1 - y) / n_sim)
low_perc = np.concatenate((scipy.stats.norm.ppf(0.005, y[:-1], std[:-1]), (1.0, )))
high_perc = np.concatenate((scipy.stats.norm.ppf(0.995, y[:-1], std[:-1]), (1.0, )))
ax[1].fill_between(L, low_perc, high_perc, color='slategray', alpha=0.5, label='__nolegend__')

# ##############################################################################
#  RANK DISTRIBUTIONS
# ##############################################################################
for g, d in data.groupby('param'):
    if g == 'fc_sigma':
        z = 100
    else:
        z = 101
    x, y = np.sort(d['rank']), np.arange(0, len(d), 1) / len(d)
    ax[1].step(x, y, lw=2, color=param_colors[g], label=legend[g], zorder=z)
leg = ax[1].legend(title='parameter', fontsize=8)
leg.get_title().set_fontsize(8)
plt.tight_layout()
plt.savefig('../../figures/Chure2019_FigS8_empirical_F_sbc_sensitivity.pdf', 
            bbox_inches='tight')

