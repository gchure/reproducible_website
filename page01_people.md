---
layout: page
title: People
description: About the authors
img: people.png 
permalink: people
sidebar: true
---


{% for author in site.data.people %}
 * {% if author.link %}[**{{author.name}}**]({{author.link}}){% else
   %}**{{author.name}}**{% endif %} {% if author.title %} \| {{author.title}} {%
   endif %} {% if author.institute %}\| {{author.institute}} {% endif %}
 {% endfor %}
