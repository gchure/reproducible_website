# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import mut.stats
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the samples
sbc_samples = pd.read_csv('../../data/Chure2019_DNA_sbc_samples.csv')

# ##############################################################################
# FIGURE 1: PRIOR DISTRIBUTIONS
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

# Add labels
for i in range(2):
    ax[i, 0].set_xlabel('DNA binding energy [$k_BT$]', fontsize=8)
    ax[i, 1].set_xlabel('$\sigma$', fontsize=8)
    ax[0, i].set_ylabel('$\propto$ probability', fontsize=8)
    ax[1, i].set_ylabel('cumulative distribution', fontsize=8)

# Add title
ax[0, 0].set_title(r'DNA binding energy $\Delta\varepsilon_{RA}$', fontsize=8,
                    y=1.08, backgroundcolor=colors['pale_yellow'])
ax[0, 1].set_title(r'standard deviation $\sigma$', fontsize=8,
                    y=1.08, backgroundcolor=colors['pale_yellow'])
axes = {'ep_RA':0, 'sigma':1}

fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# ##############################################################################
#  GROUND TRUTH DISTRIBUTIONS
# ##############################################################################
for g, d in sbc_samples.groupby(['param']):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Histogram
    hist, bins = np.histogram(d['ground_truth'], bins=15, density=True)
    hist_ax.step(bins[:-1], hist, color=colors['blue'], lw=1, label='ground truth')
    hist_ax.fill_between(bins[:-1], hist, step='pre', color=colors['blue'], alpha=0.4)
    
    # ECDF
    x, y = np.sort(d['ground_truth']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['blue'], lw=1, label='ground truth')

# ##############################################################################
#  SBC DISTRIBUTIONS
# ##############################################################################
for g, d in sbc_samples.groupby(['param']):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Histogram
    hist, bins = np.histogram(d['post_mean'], bins=15, density=True)
    hist_ax.step(bins[:-1], hist, color=colors['red'], lw=1, label='inferred')
    hist_ax.fill_between(bins[:-1], hist, step='pre', color=colors['red'], alpha=0.4)
    
    # ECDF
    x, y = np.sort(d['post_mean']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['red'], lw=1, label='inferred')

# ##############################################################################
# LEGENDS
# ##############################################################################
ax[0, 0].legend(fontsize=8, handlelength=0.5)
plt.tight_layout()
plt.savefig('../../figures/FigS3_DNA_prior_recovery.pdf', bbox_inches='tight')

# ##############################################################################
# FIGURE 2: SENSITIVITY
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(4, 2.5), sharex=True, sharey=True)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

# Add labels
for i in range(2):
    ax[i].set_xlabel('posterior shrinkage $s$', fontsize=8)
ax[0].set_ylabel('posterior $z$-score', fontsize=8)

# Add titles
ax[0].set_title(r'DNA binding energy $\Delta\varepsilon_{RA}$', fontsize=8, y=1.08,
                backgroundcolor=colors['pale_yellow'])
ax[1].set_title('standard deviation $\sigma$', fontsize=8, 
                backgroundcolor=colors['pale_yellow'], y=1.08)

# Add panel labels
fig.text(0.02, 0.9, '(A)', fontsize=8)
fig.text(0.53, 0.9, '(B)', fontsize=8)

# Adjust scaling
for i in range(2):
    ax[i].set_xlim([0, 1.05])
    ax[i].set_ylim([-5, 5])

# Assign axes
axes = {'ep_RA':0, 'sigma':1}
# ##############################################################################
# SBC DATA
# ##############################################################################
for g, d in sbc_samples.groupby(['param']):
    _ax = ax[axes[g]]
    _ax.plot(d['shrinkage'], d['z_score'], '.', color=colors['red'], ms=1)

plt.tight_layout()
plt.savefig('../../figures/Chure2019_FigS4_DNA_sbc_sensitivity.pdf')

# ##############################################################################
# FIGURE 3: RANK DISTRIBUTION
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('rank statistic', fontsize=8)


# Formatting and scale
for i in range(2):
    ax[0, i].set_xlim([0, 800])
    ax[1, i].set_ylabel('cumulative distribution', fontsize=8)
    ax[0, i].set_ylabel('counts', fontsize=8)

# Add labels
ax[0,0].set_title(r'DNA binding energy $\Delta\varepsilon_{RA}$', y=1.08,
                backgroundcolor=colors['pale_yellow'], fontsize=8)
ax[0,1].set_title(r'standard deviation $\sigma$', y=1.08,
                backgroundcolor=colors['pale_yellow'], fontsize=8)
fig.text(0.0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# Define the axes
axes = {'ep_RA': 0, 'sigma':1}
# ##############################################################################
# TRUE UNIFORM DISTRIBUTION
# ##############################################################################
n_sim = sbc_samples.sim_idx.max()
bins = 20
L = np.arange(0, n_sim, 1)
R = 800

# Bounds for histogram 99% 
low = scipy.stats.binom.ppf(0.005, R, 2 * bins/R)
high = scipy.stats.binom.ppf(0.995, R, 2* bins/R)

# Envelope of cdf 99%
y = scipy.stats.randint.cdf(L, 0, R)
std = np.sqrt(y * (1 - y) / n_sim)
low_perc = np.concatenate((scipy.stats.norm.ppf(0.005, y[:-1], std[:-1]), (1.0, )))
high_perc = np.concatenate((scipy.stats.norm.ppf(0.995, y[:-1], std[:-1]), (1.0, )))

# ##############################################################################
#  DATA DISTRIBUTION
# ##############################################################################
for g, d in sbc_samples.groupby('param'):
    hist_ax = ax[0, axes[g]]
    ecdf_ax = ax[1, axes[g]]

    # Bin the histogram
    _ = hist_ax.hist(d['rank'], bins=bins, color=colors['red'], 
                edgecolor='k')

    # Percentile bounds
    _ = hist_ax.fill_between(L, low, high, color='slategray', alpha=0.4, zorder=100)

    # ECDF
    x, y = np.sort(d['rank']), np.arange(0, len(d), 1) / len(d)
    ecdf_ax.step(x, y, color=colors['red'])

    # Percentile_bounds
    ecdf_ax.fill_between(L, low_perc, high_perc, color='slategray', alpha=0.4)

plt.tight_layout()
plt.savefig('../../figures/Chure2019_FigS5_DNA_sbc_rank_distribution.pdf', bbox_inches='tight')