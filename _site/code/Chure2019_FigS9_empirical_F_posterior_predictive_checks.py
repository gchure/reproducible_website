# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mut.bayes
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/Chure2019_compiled_data.csv', comment='#')
data = data[(data['mutant']=='Y20I-Q294V') & (data['IPTGuM']==50)]
model = mut.bayes.StanModel('../stan/Chure2019_empirical_F_inference.stan')
_, samples = model.sample(dict(N=len(data), foldchange=data['fold_change']))
samples.sort_values(by='lp__', inplace=True)
v_cmap = sns.color_palette('viridis', n_colors=len(samples))
samples['color'] = v_cmap

dfs = []
mu = samples['fc_mu'].values
sig = samples['fc_sigma'].values
for i in range(len(samples)):
    n_rand = np.random.normal(mu[i], sig[i], len(data))
    df = pd.DataFrame([])
    df['fold_change'] = n_rand
    df['draw'] = i
    df['mutant'] = 'Y20I-Q294V'
    df['IPTGuM'] = 50
    dfs.append(df)
df = pd.concat(dfs)

percs = [99, 95, 80, 50, 20, 10, 5]
perc_cmap = sns.color_palette('viridis', n_colors=(len(percs)))
zorder = [11, 12, 13, 14, 15, 16, 17]
z = {p:z for p, z in zip(percs, zorder)}
c = {p:c for p, c in zip(percs, perc_cmap)}

# Group by the ECDFS
dfs = []
for g, d in df.groupby(['draw']):
    d = d.copy()
    d = d.sort_values('fold_change')
    d['y'] = np.arange(0, len(d), 1) / len(d)
    dfs.append(d)
cdf_data = pd.concat(dfs)

# And compute the percentiles
perc_df = pd.DataFrame([])
for g, d in cdf_data.groupby('y'):
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

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(6, 2.5))
gs = gridspec.GridSpec(2, 5)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])
ax4 = fig.add_subplot(gs[:, 3:])
ax = [ax1, ax2, ax3, ax4]
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# Disable axes where needed
ax1.set_xticks([])
ax1.set_yticks([])
ax3.set_yticks([])

# Add labels
ax2.set_xlabel('$\mu$', fontsize=8)
ax2.set_ylabel('$\sigma$', fontsize=8)
ax3.set_xlabel('$\sigma$', fontsize=8)
ax4.set_xlabel('fold-change', fontsize=8)
ax4.set_ylabel('cumulative distribution', fontsize=8)

# Adjust limits
ax4.set_xlim([0, 1])

# Add panel labels
fig.text(0.05, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# ##############################################################################
# JOINT DISTRIBUTION
# ##############################################################################
ax2.scatter(samples['fc_mu'], samples['fc_sigma'], c=samples['color'], marker='.', s=1)

# ##############################################################################
#  MARGINAL DISTRIBUTIONS
# ##############################################################################
ax1.hist(samples['fc_mu'], bins=20, histtype='stepfilled', color=v_cmap[10],
        alpha=0.5)
ax3.hist(samples['fc_sigma'], bins=20, histtype='stepfilled', color=v_cmap[10],
        alpha=0.5)

# ##############################################################################
# PERCENTILE
# ##############################################################################
for g, d in df.groupby('percentile'):
   ax4.fill_betweenx(np.linspace(0, 1, len(d)), d['fc_low'], d['fc_high'], 
                    color=c[g], zorder=z[g], alpha=0.5, label=int(g))

# ##############################################################################
# DATA
# ##############################################################################    
x = np.sort(data['fold_change'])
y = np.arange(0, len(x), 1) / len(x)
ax4.step(x, y, 'k-', markerfacecolor='w', ms=3, label='__nolegend__', zorder=999)
ax4.plot(x, y, 'ko', markerfacecolor='w', ms=3, label='data',zorder=1000)
leg = ax4.legend(title='percentile', fontsize=8, loc='upper left', handlelength=1)
leg.get_title().set_fontsize(8)
plt.savefig('../../figures/FigS9_empirical_F_posterior_predictive.pdf', 
            bbox_inches='tight')



