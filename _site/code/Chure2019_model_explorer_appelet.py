# -*- coding: utf-8 -*-
import numpy as np
from bokeh.themes import Theme
import mut.viz
import mut.thermo
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
from bokeh.embed import components
pboc = mut.viz.color_selector('pboc')
bokeh.plotting.output_file("../../figures/model_explorer.html")

# Define parameter ranges
c_range = np.logspace(-8, 8, 500)
subsamp = list(np.arange(0, len(c_range), 25))
bohr_range = np.linspace(-20, 20, 500)

# Set the reference induction profile
ref_arch = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, effector_conc=c_range,
                                    ka=139, ki=0.53, ep_ai=4.5, 
                                    n_sites=2)
ref_fc = ref_arch.fold_change()
ref_bohr = ref_arch.bohr_parameter()
ref_delta_bohr = ref_bohr - ref_bohr
# Define the source
source = ColumnDataSource(data=dict(c=c_range, ref_fc=ref_fc, 
                                        mut_fc=ref_fc, mut_bohr=ref_bohr,
                                        mut_delta_bohr=ref_delta_bohr, ref_bohr=ref_bohr))
view = CDSView(source=source, filters=[IndexFilter(subsamp)])


# Instantiate the figure canvas
p_fc = bokeh.plotting.figure(width=400, height=300, x_axis_type='log', 
                            x_axis_label='IPTG [µM]', y_axis_label='fold-change',
                            x_range=[1E-2, 1E5], y_range=[-0.1, 1.1],
                            title='    INDUCTION PROFILE    ')
p_bohr = bokeh.plotting.figure(width=425, height=300,
                                x_axis_label='free energy [kT]', 
                                y_axis_label='fold-change',
                                x_range=[-15, 15], y_range=[-0.1, 1.1],
                                title='    PHENOTYPIC DATA COLLAPSE    ')
p_delBohr = bokeh.plotting.figure(width=425, height=300,
                                  x_axis_label='IPTG [µM]', y_axis_label='∆F [kT]',
                                  x_axis_type='log',
                                  y_range=[-15, 15],
                                  title='    CHANGE IN FREE ENERGY    ')

# Plot fold-change values
p_fc.line(x='c', y='ref_fc', source=source, color=pboc['blue'], line_width=2,
          legend='reference strain')
p_fc.line(x='c', y='mut_fc', source=source, color=pboc['red'], line_width=2,
          legend='mutant strain')
p_fc.circle(x='c', y='mut_fc', source=source, view=view,
            color=pboc['red'], line_width=2, size=8, legend='mutant strain')


# Plot the data collapse
p_bohr.line(bohr_range, (1 + np.exp(-bohr_range))**-1, color='black', 
            line_width=2, legend='master curve')
p_bohr.circle(x='mut_bohr', y='mut_fc', source=source, view=view,
            color=pboc['red'], size=8, legend='mutant strain')

# Pot the delta bohr
dbohr_ref = bokeh.models.Span(location=0, dimension='width', line_color='black',
                               line_width=2)
p_delBohr.line(x='c', y='mut_delta_bohr', source=source, color=pboc['red'], 
            legend='mutant strain', line_width=2)
p_delBohr.circle(x='c' , y='mut_delta_bohr', view=view, source=source,
                color=pboc['red'],  size=8, legend='mutant strain')
dbohr_ref = bokeh.models.Span(location=0, dimension='width', line_color='black',
                               line_width=2)
p_delBohr.renderers.extend([dbohr_ref])


# Position legends
p_bohr.legend.location = 'top_left'
p_delBohr.legend.location = 'top_left'
p_fc.legend.location = 'top_left'


# #######################
# CONTROLS
# #######################
# Reference controls
ref_epRA_slider = Slider(start=-20, end=-2, value=-13.9, step=0.1,
title='DNA binding energy [kT]', bar_color=pboc['blue'])
ref_R_slider = Slider(start=1, end=2000, value=260, step=1,
title='Repressors per cell', bar_color=pboc['blue'])
ref_epAI_slider = Slider(start=-10, end=10, value=4.5, step=0.1,
title='Allosteric state energy difference [kT]', bar_color=pboc['blue'])
ref_ka_slider = Slider(start=-2, end=3.7, value=np.log10(139), step=0.1,
    title='Log10(Ka/1µM)', bar_color=pboc['blue'])
ref_ki_slider = Slider(start=-3, end=2.7, value=np.log10(0.53), step=0.01,
    title='Log10(Ki/1µM)', bar_color=pboc['blue'])
ref_n_slider = Slider(start=1, end=10, value=2, step=1,
    title='Number of allosteric sites', bar_color=pboc['blue'])

# Mutant controls
mut_epRA_slider = Slider(start=-20, end=-2, value=-13.9, step=0.1,
title='DNA binding energy [kT]', bar_color=pboc['light_red'])
mut_R_slider = Slider(start=1, end=2000, value=260, step=1,
title='Repressors per cell', bar_color=pboc['light_red'])
mut_ka_slider = Slider(start=-2, end=3.7, value=np.log10(139), step=0.01,
    title='Log10(Ka / 1µM)', bar_color=pboc['light_red'])
mut_ki_slider = Slider(start=-3, end=2.7, value=np.log10(0.53), step=0.01,
    title='Log10(Ki / 1µM)', bar_color=pboc['light_red'])
mut_epAI_slider = Slider(start=-10, end=10, value=4.5, step=0.1,
    title='Allosteric state energy difference [kT]', bar_color=pboc['light_red'])
mut_n_slider = Slider(start=1, end=10, value=2, step=1,
    title='Number of allosteric sites', bar_color=pboc['light_red'])
    
callback_args = {'source':source,
                 'refepRA':ref_epRA_slider,
                 'refR': ref_R_slider,
                 'refKa': ref_ka_slider,
                 'refKi': ref_ki_slider,
                 'refepAI':ref_epAI_slider,
                 'ref_n_sites':ref_n_slider,
                 'mutepRA':mut_epRA_slider,
                 'mutR': mut_R_slider,
                 'mutKa': mut_ka_slider,
                 'mutKi': mut_ki_slider,
                 'mutepAI':mut_epAI_slider,
                 'mut_n_sites':mut_n_slider}

reset_ref = CustomJS(args=callback_args, code="""
        refepRA.value = -13.9;
        refR.value = 260;
        refKa.value = Math.log10(139);
        refKi.value = Math.log10(0.53);
        refepAI.value = 4.5;
        ref_n_sites.value = 2;
        """)
reset_mut = CustomJS(args=callback_args, code="""
        mutepRA.value = refepRA.value;
        mutR.value = refR.value;
        mutKa.value = refKa.value;
        mutKi.value = refKi.value;
        mutepAI.value = refepAI.value;
        mut_n_sites.value = ref_n_sites.value;
        """)

callback  = CustomJS(args=callback_args, code="""
                // Define constants
                var data = source.data;
                var c = data['c'];
                // Reference parameters
                var ref_fc = data['ref_fc'];
                var ref_ep_r = refepRA.value;
                var ref_R = refR.value;
                var ref_ka = Math.pow(10,refKa.value);
                var ref_ki = Math.pow(10, refKi.value);
                var ref_ep_ai= refepAI.value;
                var ref_n = ref_n_sites.value;

                // Mutant parameters
                var mut_fc = data['mut_fc'];
                var mut_bohr = data['mut_bohr'];
                var mut_delta_bohr = data['mut_delta_bohr'];
                var mut_ep_r = mutepRA.value;
                var mut_R = mutR.value;
                var mut_ka = Math.pow(10, mutKa.value);
                var mut_ki = Math.pow(10, mutKi.value);
                var mut_ep_ai= mutepAI.value;
                var mut_n = mut_n_sites.value;

                // Define functions for calculating various quantities
                function computePact(c, ka, ki, ep_ai, n) {
                    var numer = Math.pow(1 + c / ki, n);
                    var denom = Math.pow(1 + c / ka, n);
                    return Math.pow(1 + Math.exp(-ep_ai) * numer / denom, -1);
                }

                function computeFoldChange(R, ep_r, c, ka, ki, ep_ai, n) {
                    var pact = computePact(c, ka, ki, ep_ai, n);
                    // Note that 4.6E6 is hard coded as the number of nonspecific 
                    // sites
                    return Math.pow(1 + pact * (R / 4600000) * Math.exp(-ep_r), -1);
                }

                function computeBohrParameter(R, ep_r, c, ka, ki, ep_ai, n) {
                    var pact = computePact(c, ka, ki, ep_ai, n);
                    // Note that 4.6E6 is hard coded as the number of nonspecific 
                    // sites
                    return -Math.log(pact) - Math.log(R/4600000) + ep_r;
                }

                // Evaluate the fold-change of the reference and perturbed state
                for (var i = 0; i < c.length; i++) {
                    ref_fc[i] = computeFoldChange(ref_R, ref_ep_r, c[i], ref_ka, 
                                                  ref_ki, ref_ep_ai, ref_n);
                    mut_fc[i] = computeFoldChange(mut_R, mut_ep_r, c[i], mut_ka, 
                                                  mut_ki, mut_ep_ai, mut_n);
                    mut_bohr[i] =  computeBohrParameter(mut_R, mut_ep_r, c[i],
                                                        mut_ka, mut_ki, mut_ep_ai,
                                                        mut_n);
                    var ref_bohr = computeBohrParameter(ref_R, ref_ep_r, c[i],
                                                        ref_ka, ref_ki, ref_ep_ai,
                                                        ref_n);
                    source.data['ref_bohr'][i] = ref_bohr;
                    mut_delta_bohr[i] = mut_bohr[i] - ref_bohr;
                }
                source.change.emit();
                """)

# Define the buttons
ref_reset = Button(label='double-click to reset wild type', callback=reset_ref)
mut_reset = Button(label='double-click to set mutant to wild type', callback=reset_mut)
ref_reset.js_on_click(callback)
mut_reset.js_on_click(callback)

# Assemble controls
ref_controls = [ref_reset, ref_R_slider,  ref_epRA_slider, ref_ka_slider, ref_ki_slider,
                ref_epAI_slider, ref_n_slider]
mut_controls = [mut_reset, mut_R_slider,  mut_epRA_slider, mut_ka_slider, mut_ki_slider,
                mut_epAI_slider, mut_n_slider]

for rc, mc in zip(ref_controls[1:], mut_controls[1:]):
    rc.callback = callback
    mc.callback = callback
ref_inputs = widgetbox(ref_controls)
mut_inputs = widgetbox(mut_controls)
layout = bokeh.layouts.layout([[ref_inputs, mut_inputs, p_fc], 
                               [p_bohr, p_delBohr]])

# #############################
# THEME DETAILS
# ############################
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
bokeh.io.save(layout)

