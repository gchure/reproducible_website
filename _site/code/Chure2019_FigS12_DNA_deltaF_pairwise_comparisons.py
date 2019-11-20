# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the empirical F data
data = pd.read_csv('../../data/Chure2019_empirical_F_statistics.csv')
data = data[data['class']=='DNA']
epRA_stats = pd.read_csv('../../data/Chure2019_DNA_binding_energy_summary.csv')



# ##############################################################################
#  FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(4, 3, figsize=(7, 4), sharex=True, sharey=True)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xscale('symlog', linthreshx=2E-2)
    a.set_ylim([-8, 8])

for i in range(3):
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=6)

for i in range(4):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=6)

# Define the axes. 
cols = {'Q21A':1, 'Q21M':0, 'Y20I':2}
rows = {60:0, 124:1, 260:2, 1220:3}

# Add titles and row labels. 
for m, a in cols.items():
    ax[0, a].set_title(m, fontsize=6, y=1.08, backgroundcolor=colors['pale_yellow'])
    
for r, a in rows.items():
    ax[a, 0].text(-0.3, 0.57, f'R = {int(r)}', fontsize=6, rotation='vertical',
                backgroundcolor=colors['pale_yellow'],
                transform=ax[a, 0].transAxes)

# Define the repressor colors
rep_colors = {60:colors['purple'], 124:colors['blue'], 
              260:colors['red'], 1220:colors['green']}

# ##############################################################################
# PLOT LEGEND ENTRIES
# ##############################################################################
for r, c in rep_colors.items():
    ax[0, 0].plot([], [], '.', ms=5, color=c, label=int(r))
leg = ax[0,  0].legend(title='repressors per cell', fontsize=6, ncol=4,
                    columnspacing=0.2)
leg.get_title().set_fontsize(6)

# ##############################################################################
# PLOT DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'repressors', 'IPTGuM']):
    dbohr = d[d['parameter']=='delta_bohr']
    fc_mu = d[d['parameter']=='fc_mu']['median'].values[0]
    fc_sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc_mu > fc_sig) & (1 - fc_mu > fc_sig):
        for r, a in rows.items():
            _ax = ax[a, cols[g[0]]]
            if g[1] == r:
                face = 'w'
                zorder=1000
            else:
                face = rep_colors[g[1]]
                zorder=100
            _ax.plot(dbohr['IPTGuM'], dbohr['median'], '.', color=rep_colors[g[1]],
            markerfacecolor=face, markeredgewidth=0.75, label='__nolegend__',
            zorder=zorder)
            _ax.vlines(dbohr['IPTGuM'], dbohr['hpd_min'], dbohr['hpd_max'], lw=0.75,
                color=rep_colors[g[1]], label='__nolegend__', zorder=1)


# ##############################################################################
#  PREDICTIONS
# ##############################################################################
for c, ca in cols.items():
    for r, ra in rows.items():
        _ax = ax[ra, ca]
        epRA = epRA_stats[(epRA_stats['parameter']=='ep_RA') & 
        (epRA_stats['mutant']==c) & (epRA_stats['repressors']==r)]
        
        _ax.fill_between([0, 1E4], epRA['hpd_min'] - constants['O2'], 
                        epRA['hpd_max'] - constants['O2'], alpha=0.75, 
                        color=rep_colors[r])
plt.subplots_adjust(hspace=0.05, wspace=0.03)
plt.savefig('../../figures/Chure2019_FigS12_DNA_deltaF_pairwise_comparison.pdf', 
            bbox_inches='tight')
