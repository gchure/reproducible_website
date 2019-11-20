# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.bayes
import mut.thermo
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the data and sampling statistics
data = pd.read_csv('../../data/Chure2019_empirical_F_simulated_data.csv')
razo_data = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
razo_delF = pd.read_csv('../../data/RazoMejia2018_empirical_F_statistics.csv')
stats = pd.read_csv('../../data/Chure2019_empirical_F_simulated_data_statistics.csv')
stats['draw'] = stats.groupby(['true_bohr']).ngroup()
stats
# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
ax[0, 0].set_xlabel('free energy [$k_BT$]', fontsize=6)
ax[1, 0].set_xlabel('free energy [$k_BT$]', fontsize=6)
ax[0, 0].set_ylabel('fold-change', fontsize=6)
ax[1, 0].set_ylabel('fold-change', fontsize=6)
ax[0, 1].set_xlabel('true free energy [$k_BT$]', fontsize=6)
ax[1, 1].set_xlabel('predicted free energy [$k_BT$]', fontsize=6)
ax[1, 1].set_ylabel('inferred free energy [$k_BT$]', fontsize=6)
ax[0, 1].set_ylabel('inferred free energy [$k_BT$]', fontsize=6)
ax[0, 0].set_ylim([-0.25, 1.45])
ax[1, 0].set_ylim([-0.25, 1.45])
ax[1, 1].set_xlim([-10, 10])
ax[1, 0].set_xlim([-10, 10])

# Panel labels. 
fig.text(0.0, 1, '(A)', fontsize=8)
fig.text(0.5, 1, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)
# ##############################################################################
# AGREEMENT AND MASTER CURVE
# ##############################################################################
bohr_range = np.linspace(-15, 15, 200)
for i in range(2):
    ax[i, 0].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', 
                        label='__nolegend__')
    ax[i, 1].plot(bohr_range, bohr_range, 'k-', label='perfect agreement')

# ##############################################################################
# SIMULATED DATA
# ##############################################################################
for g, d in data.groupby(['draw']):
    if g%10 == 0:
        if g == 0:
            label = 'true fold-change'
        else:
            label = '__nolegend__'
        ax[0, 0].plot(d['bohr'], d['fc_mu'], '.', ms=5, 
                         color='k', label=label, zorder=100,
                         markerfacecolor='w', markeredgewidth=0.75)
        if g == 0:
            label = 'simulated data'
        else:
            label = '__nolegend__'
        ax[0, 0].plot(d['bohr'], d['fold_change'], '.', ms=3, 
                        color=colors['red'], zorder=99,
                        label=label)

# ##############################################################################
# INFERRED DATA
# ##############################################################################
labels = []
for g, d in stats.groupby('draw'):
    _fc = d[d['parameter']=='fc_mu']
    _sig = d[d['parameter']=='fc_sigma']
    _F = d[d['parameter']=='empirical_bohr']
    if g%10 == 0:
        if g == 0:
            legend = 'inferred fold-change'
        else:
            legend = '__nolegend__'
        ax[0, 0].plot(d['true_bohr'].values[0], _fc['median'], '.', ms=5, 
                        color=colors['blue'], label=legend, zorder=1001,
                        markerfacecolor=colors['blue'],
                        markeredgewidth=0.75)
        ax[0, 0].vlines(d['true_bohr'].values[0], _fc['hpd_min'], _fc['hpd_max'], 
                    lw=0.75, color=colors['blue'], label='__nolegend__')

        if (_fc['median'].values[0] < _sig['median'].values[0]):
            _color = colors['purple']
            face = colors['purple']
            label = 'µ < σ'
        elif (1 - _fc['median'].values[0] < _sig['median'].values[0]):
            _color = colors['green']
            face = colors['green']
            label = '1 - µ < σ'
        else:
            _color = colors['blue']
            face = colors['blue']
            label = 'µ > σ; 1 - μ > σ'
        if label not in labels:
            labels.append(label)
        else:
            label = '__nolegend__'

            ax[0, 1].plot(d['true_bohr'].values[0], _F['median'], '.', 
                        color=_color, label=label, ms=5, markerfacecolor=face, markeredgewidth=0.75)
            ax[0,1].vlines(d['true_bohr'].values[0], _F['hpd_min'], 
                        _F['hpd_max'], lw=0.75, color=_color,     label='__nolegend__') 

# ##############################################################################
# INDUCTION PAPER COLLAPSE
# ##############################################################################
bohr = mut.thermo.SimpleRepression(R=razo_data['repressors'] * 2, 
                                   ep_r=razo_data['binding_energy'],
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   ep_ai=constants['ep_AI'],
                                   effector_conc=razo_data['IPTG_uM']).bohr_parameter()
razo_data['bohr'] = bohr
ax[1, 0].plot(razo_data['bohr'], razo_data['fold_change_A'], '.', color=colors['red'],
            ms=1, label='Razo-Mejia et al. 2018', zorder=1, alpha=0.75)

# ##############################################################################
# INFERRED INDUCTION DATA
# ##############################################################################
labels = []
for g, d in razo_delF.groupby(['pred_bohr']):
    fc_mu = d[d['parameter']=='fc_mu']
    fc_sig = d[d['parameter']=='fc_sigma']
    delF = d[d['parameter']=='empirical_bohr']
    if (fc_mu['median'].values[0] < fc_sig['median'].values[0]):
        color=colors['purple']
        label = 'µ < σ'
    elif (1 - fc_mu['median'].values[0] < fc_sig['median'].values[0]):
        color = colors['green']
        label = '1 - µ < σ'
    else:
        color = colors['blue']
        label =  'µ > σ ; 1 - µ > σ'
    if label not in labels:
        labels.append(label)
    else:
        label = '__nolegend__'
    ax[1, 0].plot(fc_mu['pred_bohr'], fc_mu['median'], '.', ms=3, color=color,
            label=label, 
            markeredgewidth=0.75)
    ax[1, 0].vlines(fc_mu['pred_bohr'], fc_mu['hpd_min'], fc_mu['hpd_max'], 
                lw=0.5, color=color, label='__nolegend__')
        
    ax[1, 1].plot(delF['pred_bohr'], delF['median'], '.', color=color, label=label,
        ms=2)
    ax[1, 1].vlines(delF['pred_bohr'], delF['hpd_min'], delF['hpd_max'], lw=0.5,
    color=color, label='__nolegend__')

# ##############################################################################
# LEGEND AND SAVING
# ##############################################################################
ax[0, 0].legend(loc='upper left', fontsize=6)
ax[0, 1].legend(loc='upper left', fontsize=6)
ax[1, 0].legend(loc='upper left', fontsize=6)
ax[1, 1].legend(loc='upper left', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/Chure2019_FigS10_empirical_F_systematic_error.pdf', 
            bbox_inches='tight')
