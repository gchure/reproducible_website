# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot  as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import mut.viz
import mut.thermo
import mut.stats
import scipy.stats
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
ppc_data = pd.read_csv('../../data/Chure2019_IND_prior_predictive_checks.csv')

ep_a_unique = ppc_data['ka'].unique()
ep_i_unique = ppc_data['ki'].unique()
ep_ai_unique = ppc_data[ppc_data['model']=='KaKi_epAI']['ep_ai'].unique()
sig_unique = ppc_data['sigma'].unique()

# ##############################################################################
# FIGURE INSTANTIATION #
# ##############################################################################
fig = plt.figure(figsize=(7, 4), dpi=120)
gs = gridspec.GridSpec(4, 9)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[0:2,2:4])
ax3 = fig.add_subplot(gs[2:4, 0:2])
ax4 = fig.add_subplot(gs[2:4, 2:4])
ax5 = fig.add_subplot(gs[0:2, 6:])
ax6 = fig.add_subplot(gs[2:4, 6:])
ax = [ax1, ax2, ax3, ax4, ax5, ax6]
for a in ax:
   a.xaxis.set_tick_params(labelsize=6)
   a.yaxis.set_tick_params(labelsize=6)

for a in ax[:4]:
   a.set_yticks([])

# Add axis labels
ax1.set_xlabel(r'$K_A$ [µM]', fontsize=8, labelpad=0.1)
ax2.set_xlabel(r'$K_I$ [µM]', fontsize=8, labelpad=0.1)
ax3.set_xlabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]', fontsize=8, labelpad=0.1)
ax4.set_xlabel(r'$\sigma$', fontsize=8, labelpad=0.1)
ax5.set_xticklabels([])
ax6.set_xlabel('IPTG [µM]', fontsize=8)
ax5.set_xlabel('IPTG [µM]', fontsize=8)
ax5.set_ylabel('fold-change', fontsize=8)
ax6.set_ylabel('fold-change', fontsize=8)

# # Set scaling
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_xlim([1E-3, 1E4])
ax2.set_xlim([1E-3, 1E4])
ax1.set_ylim([0, 0.275])
ax2.set_ylim([0, 0.275])
ax3.set_ylim([0, 0.0972])
ax3.set_xlim([-20, 20])
ax4.set_xlim([0, 1])
ax4.set_ylim([0, 2.5])
ax5.set_xscale('symlog', linthreshx=1E-2)
ax6.set_xscale('symlog', linthreshx=1E-2)
ax5.set_xlim([0, 5E3])
ax6.set_xlim([0, 5E3])

# # Set limits
models = ['$K_A$ and $K_I$ only', r'$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$']
for i, a in enumerate([ax5, ax6]):
   a.set_title(models[i], fontsize=8, y=1.04, backgroundcolor=colors['pale_yellow'])

# # Add panel labels
fig.text(0.1, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# # Define the axes
axes = {'KaKi_only':  ax5, 'KaKi_epAI': ax6}

# # ##############################################################################
# # PDF PRIOR DISTRIBUTIONS
# # ##############################################################################
k_range = np.logspace(-3, 4, 1000)
ep_range = np.linspace(-18, 18, 100)
sig_range = np.linspace(0, 1, 200)
ka_pdf = scipy.stats.norm(2, 2).pdf(np.log(k_range))
ki_pdf = scipy.stats.norm(0, 2).pdf(np.log(k_range))
epAI_pdf = scipy.stats.norm(0, 5).pdf(ep_range)
sig_pdf = scipy.stats.norm(0, 0.2).pdf(sig_range)
ax1.plot(k_range, ka_pdf, '-', color=colors['red'])
ax1.fill_between(k_range, ka_pdf, '-', color=colors['light_red'])
ax2.plot(k_range, ki_pdf, '-', color=colors['red'])
ax2.fill_between(k_range, ki_pdf, '-', color=colors['light_red'])
ax3.plot(ep_range, epAI_pdf, '-', color=colors['blue'])
ax3.fill_between(ep_range, epAI_pdf, '-', color=colors['light_blue'])
ax4.plot(sig_range, sig_pdf, '-', color=colors['red'])
ax4.fill_between(sig_range, sig_pdf, '-', color=colors['light_red'])

# ##############################################################################
# SAMPLED PRIOR DISTRIBUTIONS
# ##############################################################################
k_dist = np.random.normal(0.25, 0.008, len(ep_a_unique))
ep_dist = np.random.normal(0.09, 0.002, len(ep_a_unique))
sig_dist = np.random.normal(2.35, 0.05, len(sig_unique))
ax1.plot(ep_a_unique, k_dist, 'k.', ms=0.5, alpha=0.5)
ax2.plot(ep_i_unique, k_dist, 'k.', ms=0.5, alpha=0.5)
ax3.plot(ep_ai_unique, ep_dist, 'k.', ms=0.5, alpha=0.5)
ax4.plot(sig_unique, sig_dist, 'k.', ms=0.5, alpha=0.5 )
ax1.hlines(0.25, 1E-4, 1E4, color='w', lw=25)
ax2.hlines(0.25, 1E-4, 1E4, color='w', lw=25)
ax3.hlines(0.09, -20, 20, color='w', lw=15) 
ax4.hlines(2.35, 0, 1, color='w', lw=15)

# ##############################################################################
#  PRIOR PREDICTIVE CHECKS
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
cmap_kakionly = {p:c for p, c in zip(percs, sns.color_palette('Reds', len(percs)))}
cmap_kakiepai = {p:c for p, c, in zip(percs, sns.color_palette('Blues', len(percs)))}
zorder = {p:i for p, i in zip(percs, [10, 11, 12, 13, 14, 15, 16, 17])}

# Compute the percentiles of the simulations. 
grouped = ppc_data.groupby(['IPTGuM', 'model'])
df = pd.DataFrame([], columns=['percentile', 'IPTGuM', 'fc_low', 'fc_high', 'model'])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_draw'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g[0],
                         'fc_low':_percs[0],
                         'fc_high': _percs[1],
                         'model': g[1]},
                      ignore_index=True)

for g, d in  df.groupby(['model', 'percentile']):
   _ax = axes[g[0]]
   if g[0] == 'KaKi_only':
      cmap = cmap_kakionly
   else:
      cmap = cmap_kakiepai
   _ax.fill_between(d['IPTGuM'], d['fc_low'], d['fc_high'], color=cmap[g[1]],
                  zorder=zorder[g[1]], label = g[1])
leg = ax5.legend(title='percentile', fontsize=6, bbox_to_anchor=(-0.3, 1))
leg.get_title().set_fontsize(6)
leg = ax6.legend(title='percentile', fontsize=6, bbox_to_anchor=(-0.3, 1))
leg.get_title().set_fontsize(6)
plt.subplots_adjust(wspace=0.6, hspace=2.45)
plt.savefig('../../figures/Chure2019_FigS13_IND_prior_predictive_checks.pdf', 
         bbox_inches='tight', background='white')
