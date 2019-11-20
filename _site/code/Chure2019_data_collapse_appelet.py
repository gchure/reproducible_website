# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import bokeh.plotting
import bokeh.io
from bokeh.themes import Theme
import mut.viz
import mut.thermo
bokeh.plotting.output_file("../../figures/data_collapse.html")

constants = mut.thermo.load_constants()
constants['Oid'] = -17
pboc = mut.viz.color_selector('pboc')

# Load the data sets and tweak as necessary
old_gods = pd.read_csv('../../data/Garcia2011_Brewster2014.csv', comment='#')
old_gods['IPTGuM'] = 0
old_gods['mutant'] = 'WT'
old_gods.loc[old_gods['author']=='brewster', 'method'] = 'Microscopy'
old_gods.loc[old_gods['author']=='brewster', 'work'] = 'Brewster et al. 2014'
old_gods.loc[old_gods['author']=='garcia', 'method'] = 'Miller Assay'
old_gods.loc[old_gods['author']=='garcia', 'work'] = 'Garcia & Phillips 2011'
old_gods.rename(columns={'repressor':'repressors'}, inplace=True)
new_gods = pd.read_csv('../../data/RazoMejia_2018.csv', comment='#')
new_gods = new_gods[new_gods['repressors'] > 0]
new_gods['mutant'] = 'WT'
new_gods['method'] = 'Flow Cytometry'
new_gods['work'] = "Razo Mejia et al. 2018"
data = pd.read_csv('../../data/Chure2019_summarized_data.csv', comment='#')
data['method'] = 'Flow Cytometry'
data['work'] = "This Study"


# Correct for tetramer counts
new_gods['repressors'] *= 2
new_gods.rename(columns={'IPTG_uM':'IPTGuM', 'fold_change_A':'fold_change'},
                inplace=True)
new_gods = new_gods.groupby(['repressors', 'operator', 'IPTGuM', 'method', 'mutant', 'work'])['fold_change'].agg(('mean', 'sem')).reset_index()
# Load the statistics for each sample
epRA_stats = pd.read_csv('../../data/Chure2019_DNA_binding_energy_summary.csv')
epRA_stats = epRA_stats[epRA_stats['repressors']==260]
allo_stats = pd.read_csv('../../data/Chure2019_KaKi_epAI_summary.csv')
allo_stats = allo_stats[allo_stats['operator']=='O2']

# Define colors and glyphs
_color = {'garcia': pboc_colors['red'], 'brewster':pboc_colors['blue']}
legend = {'garcia': 'Garcia & Phillips 2011', 'brewster':'Brewster et al. 2014'}


# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
TOOLTIPS = [
    ('Free Energy [kT]', '@round_bohr'),
    ('Fold-Change', '@fc_round'),
    ('IPTG [µM]', '@IPTGuM'),
    ('Operator', '@operator'),
    ('Repressors', '@repressors'),
    ('Mutant', '@mutant'),
    ('Method', '@method'),
    ('Reference', '@work')
]
p = bokeh.plotting.figure(width=600, height=600, x_axis_label='free energy [kT]',
                          y_axis_label='fold-change', x_range=[-10, 10],
                          y_range=[-0.1 ,1.4], tooltips=TOOLTIPS)

# ##############################################################################
# COLLAPSE CURVE
# ##############################################################################
F = np.linspace(-40, 40, 300)
collapse = (1 + np.exp(-F))**-1
p.line(F, collapse, color='black', line_width=2, legend='master curve')

# ##############################################################################
# Razo-Mejia  et al. 2018
# ##############################################################################
ops = [constants[o] for o in new_gods['operator'].values]
bohr = mut.thermo.SimpleRepression(R=new_gods['repressors'], ep_r=ops,
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   ep_ai=constants['ep_AI'], 
                             effector_conc=new_gods['IPTGuM']).bohr_parameter() 
new_gods['bohr'] = bohr
new_gods['round_bohr'] = np.round(new_gods['bohr'].values, decimals=1)
new_gods['fc_round'] = np.round(new_gods['mean'].values, decimals=2)
p.circle(x='bohr', y='mean', color=pboc['dark_green'],
                 legend='Razo-Mejia et al. 2018', size=12, fill_alpha=0.5,
                 line_width=1.5, source=new_gods)

# ##############################################################################
# BREWSTER AND GARCIA DATA
# ##############################################################################
for g, d in old_gods.groupby('author'):
    d = d.copy()
    ops = [constants[o] for o in d['operator'].values]
    bohr = mut.thermo.SimpleRepression(R=d['repressors'], ep_r=ops,
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=0).bohr_parameter()
    d['bohr'] = bohr
    d['round_bohr'] = np.round(d['bohr'].values, decimals=1)
    d['fc_round'] = np.round(d['fold_change'], decimals=2)
    if g=='brewster':
        p.square(x='bohr', y='fold_change', color=_color[g], 
             legend=legend[g], size=12, fill_alpha=0.5, line_width=1.5, source=d)
    else:
        p.diamond(x='bohr', y='fold_change', color=_color[g], 
             legend=legend[g], size=12, fill_alpha=0.5, line_width=1.5, source=d)


# ##############################################################################
# DNA MUTS
# ##############################################################################
op_glyphs = {'O1':'^', 'O2':'v', 'O3':'D'}
for g, d in data[data['class']=='DNA'].groupby(['mutant']):
    d = d.copy()
    ep_RA = epRA_stats[(epRA_stats['mutant']==g) &
                (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
    bohr = mut.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=constants['Ka'], ki=constants['Ki'],
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    d['bohr']  = bohr
    d['round_bohr'] = np.round(d['bohr'].values, decimals=1)
    d['fc_round'] = np.round(d['mean'].values, decimals=2)
    p.triangle(x='bohr', y='mean', color=mut_colors[g], legend=g, size=12, fill_alpha=0.5,
    fill_color='white', line_width=1.5, source=d)

# ##############################################################################
#  IND MUTS
# ##############################################################################
for g, d in data[data['class']=='IND'].groupby(['mutant']):
    d = d.copy()
    _stats = allo_stats[(allo_stats['mutant']==g)]
    ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
    ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    ops = [constants[o] for o in d['operator'].values]
    bohr = mut.thermo.SimpleRepression(R=d['repressors'], ep_r=ops, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    d['bohr'] = bohr
    d['round_bohr'] = np.round(d['bohr'].values, decimals=1)
    d['fc_round']  = np.round(d['mean'], decimals=2)
    p.hex(x='bohr', y='mean', color=mut_colors[g], legend=g, size=12, fill_alpha=0.5,
    fill_color='white', line_width=1.5, source=d)


# ##############################################################################
# DBL MUTS
# ##############################################################################
for g, d in data[data['class']=='DBL'].groupby(['mutant']):
    d = d.copy()
    ep_RA = epRA_stats[(epRA_stats['mutant']==g.split('-')[0]) & 
                           (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
    _stats = allo_stats[(allo_stats['mutant']==g.split('-')[1])]
    ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
    ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    bohr = mut.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
    d['bohr'] = bohr
    d['round_bohr'] = np.round(bohr, decimals=1)
    d['fc_round'] = np.round(d['mean'].values, decimals=2)
    p.circle_x(x='bohr', y='mean', color=mut_colors[g],
                legend=g, size=12, fill_alpha=0.5, fill_color='white', 
                line_width=1.5, source=d)


# ##############################################################################
# THEME DETAILS
# ##############################################################################
p.legend.location = 'top_left'
p.legend.click_policy = 'hide'
p.legend.label_text_font_size = '8pt'
p.legend.spacing = 1 
p.legend.padding = 3
p.legend.label_standoff = 2
theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#E3DCD0',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'axis_line_color': "white",
            'major_tick_in': 7,
            'major_tick_line_width': 2.5,
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#E3DCD0',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.show(p)
bokeh.io.save(p)