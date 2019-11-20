---
layout: page
title: Code and Data
img: code.png # Add image post (optional)
permalink: code
sidebar: true
---

---

## Executing The Code
This page should provide information of the computational environment used in
your work along with links to all of the code

## Computational Environment
All analysis and data processing was performed with the following software
configurations.

```
# Python Version
CPython 3.6.7
IPython 7.1.1

# Package Versions
bokeh==1.0.4
fcsparser==0.2.0
numpy==1.14.2
matplotlib==3.0.1
scipy==1.1.0
seaborn==0.9.0
pandas==0.23.4
tqdm==4.28.1
pystan==2.18.0.0
python_frontmatter==0.4.5
PyYAML==5.1

# System Information
compiler   : GCC 4.2.1 Compatible Clng 4.0.1 (tagss/RELEASE_401/final)
system     : Darwin
release    : 18.2.0
machine    : x86_64
processor  : i386
CPU cores  : 4
interpreter: 64bit
```

## Description of specific software module
Pug gluten-free quinoa crucifix, you probably haven't heard of them activated
charcoal meditation viral kinfolk. Schlitz mlkshk gochujang four loko
drinking vinegar farm-to-table. Edison bulb tumblr venmo post-ironic, raw
denim retro cred. Tattooed literally copper mug ethical enamel pin vinyl.
Trust fund pickled skateboard, organic poutine copper mug gentrify keytar
authentic cloud bread 8-bit. Stumptown ugh bicycle rights, cliche cronut
freegan dreamcatcher keffiyeh drinking vinegar single-origin coffee bushwick.

## Figure Generation
This section should link to all code/data used in the main text figures.

{% for fig in site.data.main_figs %}
<article class="post">

<a class="post-thumbnail" style="background-image: url({{site.baseurl}}/assets/img/{{fig.pic}})" href="{{site.baseurl}}/figures/{{fig.pdf}}"> </a>

<div class="post-content">
<b class="post-title"><a href="{{site.baseurl}}/code/{{fig.file}}">{{fig.title}}</a></b>
<p> {{fig.desc}}</p>

<i>Necessary Data Sets </i><br/>
{% for ds in fig.req %}
<a style="font-size: 0.9em;" href="{{site.baseurl}}/data/{{ds.dataset}}"> - {{ds.title}} </a><br/>
{% endfor %}
</div>
</article>
{%endfor%}