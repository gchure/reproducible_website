---
layout: page
title: Code and Data
img: code.png # Add image post (optional)
permalink: code
sidebar: true
---

---

## Description
{{site.data.code.description}}

{% if site.data.code.environment %}
## Computational Environment
```
{{site.data.code.environment}}
```

{% endif %}

## Figure Generation

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