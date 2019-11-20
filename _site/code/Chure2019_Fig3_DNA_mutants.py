# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load in the data sets
data = pd.read_csv('../../data/Chure2019_summarized_data.csv', comment='#')
data = data[data['class']=='DNA']
stats = pd.read_csv('../../data/Chure2019_DNA_binding_energy_summary.csv')

stats = stats[stats['repressors']==260]
bohr = pd.read_csv('../../data/Chure2019_empirical_F_statistics.csv')
empirical_bohr = bohr[bohr['class']=='DNA']

# Define some plotting constants. 
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
bohr_range = np.linspace(-8, 8, 200)
F = (1 + np.exp(-bohr_range))**-1

rep_colors = {60:colors['purple'], 124:colors['blue'], 
              260:colors['red'], 1220:colors['dark_green']}


# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(3, 3, figsize=(4.5, 4), dpi=150)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

for i in range(3):
    ax[0, i].set_xscale('symlog', linthreshx=0.006)
    ax[-1, i].set_xscale('symlog', linthreshx=0.006)
    ax[0, i].set_ylim([-0.2, 1.2])
    ax[1, i].set_ylim([-0.2, 1.2])
    ax[1, i].set_xlim([-8, 8])
    ax[0, i].set_xlim([-0.001, 1E4])
    ax[-1, i].set_xlim([-0.001, 1E4])
    ax[-1, i].set_ylim([-8, 8])
    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    for j in range(2):
        ax[i, j+1].set_yticklabels([])

    # Add labels   
    ax[0, i].set_xlabel("IPTG [µM]", fontsize=6)
    ax[1, i].set_xlabel("free energy [$k_BT$]", fontsize=6)
    ax[-1, i].set_xlabel("IPTG [µM]", fontsize=6)

# Add ylabels
ax[0, 0].set_ylabel('fold-change', fontsize=6)
ax[1, 0].set_ylabel('fold-change', fontsize=6)
ax[2, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=6)

# Define the axes
axes = {'Q21M':0, 'Y20I':1, 'Q21A':2}

# Add titles
for m, a in axes.items():
    ax[0, a].set_title(m, fontsize=6, backgroundcolor=colors['pale_yellow'],
    y=1.08) 

# Add panel labels
fig.text(0, 0.90, '(A)', fontsize=6)
fig.text(0, 0.63, '(B)', fontsize=6)
fig.text(0, 0.33, '(C)', fontsize=6)

# ##############################################################################
# GUIDE CURVES FOR ∆F
# ##############################################################################
for i in range(3):
    ax[-1, i].hlines(0, -0.01, 1E4, 'k', linestyle=':', lw=0.75)

# ##############################################################################
# FOLD-CHANGE CURVES
# ##############################################################################
for r, cor in rep_colors.items():
    for m, a in axes.items():
        _stats = stats[(stats['mutant']==m) & (stats['parameter']=='ep_RA')]
        _c, _ep= np.meshgrid(c_range, _stats[['hpd_min', 'hpd_max']].values)
        arch = mut.thermo.SimpleRepression(R=r, ep_r=_ep, ka=constants['Ka'],
                                           ki=constants['Ki'], 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=_c).fold_change()
        ax[0, a].fill_between(c_range, arch[0, :], arch[1, :], color=cor,
                             alpha=0.4) 

# ##############################################################################
# COLLAPSE CURVES
# ##############################################################################
for i in range(3):
    ax[1, i].plot(bohr_range, F, 'k-', lw=0.75)

# ##############################################################################
# FREE ENERGY PREDICTIONS
# ##############################################################################
for m, a in axes.items():
    _stats = stats[(stats['mutant']==m) & 
                 (stats['parameter']=='ep_RA')][['hpd_min', 'hpd_max']].values[0]
    ax[-1, a].fill_between(c_range, _stats[0] - constants['O2'], 
                         _stats[1] - constants['O2'], color=rep_colors[260], 
                    alpha=0.4)

# ##############################################################################
# FOLD-CHANGE DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'repressors']):
    if g[1] == 260:
        face = 'w'
    else:
        face = rep_colors[g[1]]
    ax[0, axes[g[0]]].errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='.',
                               color=rep_colors[int(g[1])], markerfacecolor=face,
                               ms=5, markeredgewidth=0.5, capsize=1, lw=1, linestyle='none',
                               label=int(g[1]))

# ##############################################################################
# COLLAPSE DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'repressors']):
    ep_r = stats[(stats['mutant']==g[0]) & 
            (stats['parameter']=='ep_RA')]['median'].values[0]
    bohr = mut.thermo.SimpleRepression(R=g[1], ep_r=ep_r, ka=constants['Ka'],
                                       ki=constants['Ki'], 
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=d['IPTGuM']).bohr_parameter()
    if g[1] == 260:
        face = 'w'
    else:
        face = rep_colors[g[1]]
    ax[1, axes[g[0]]].errorbar(bohr, d['mean'], d['sem'], fmt='.', 
                               linestyle='none', lw=1, capsize=1, ms=5, markeredgewidth=0.5,
                               color=rep_colors[g[1]], markerfacecolor=face)

# ##############################################################################
# INFERRED F 
# ##############################################################################
for g, d in empirical_bohr.groupby(['mutant', 'repressors', 'IPTGuM']):
    _param = d[d['parameter']=='delta_bohr']
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]

    # If the mean fold change is closer to the boundary than the sigma, do not 
    # show it (assign alpha=0)
    if (mu < sig) | (1 - mu < sig):
        color = 'slategray'
        alpha = 0
        lw = 0
        fmt = 'x'
    else:
        color = rep_colors[g[1]]
        alpha = 1 
        lw = 0.75
        fmt = '.'
    if g[1] == 260:
        face = 'w'
        zorder=1000
    elif fmt == 'x':
        zorder = 1
        face=color
    else:
        face = color
        zorder=100
    _ax = ax[-1, axes[g[0]]]
    _ax.plot(_param['IPTGuM'], _param['median'], marker=fmt, linestyle='none', 
        color=color, markerfacecolor=face, alpha=alpha, ms=5, zorder=zorder, 
        markeredgewidth=0.5)
    _ax.vlines(_param['IPTGuM'], _param['hpd_min'], _param['hpd_max'], 
            lw=lw, color=color,  alpha=alpha, zorder=zorder)

# ##############################################################################
#  LEGEND INFORMATION
# ##############################################################################
leg = ax[0, 0].legend(title='$R$', fontsize=5,handletextpad=0.1)
leg.get_title().set_fontsize(5)
plt.subplots_adjust(hspace=0.5, wspace=0.1)
plt.savefig('../../figures/Chure2019_Fig3_DNA_mutants.pdf', bbox_inches='tight')
