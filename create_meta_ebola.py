import os
from os import path
import re
import json


BASE_PATH = './Genomes/ebola'
EXT = '.fasta'


files = [f for f in os.listdir(BASE_PATH) if path.isfile(path.join(BASE_PATH, f)) and f.endswith(EXT)]
meta = []
colors = {
    'AB050936': 'aqua',
    'AF522874': 'aqua',
    'AY354458': 'blue',
    'AY729654': 'lime',
    'EU338380': 'lime',
    'FJ217161': 'silver',
    'FJ217162': 'black',
    'FJ621583': 'aqua',
    'FJ621585': 'aqua',
    'FJ968794': 'lime',
    'JN638998': 'lime',
    'JX477165': 'aqua',
    'JX477166': 'aqua',
    'KC242783': 'lime',
    'KC242784': 'fuchsia',
    'KC242785': 'fuchsia',
    'KC242786': 'fuchsia',
    'KC242787': 'fuchsia',
    'KC242788': 'fuchsia',
    'KC242789': 'fuchsia',
    'KC242790': 'fuchsia',
    'KC242791': 'purple',
    'KC242792': 'green',
    'KC242793': 'green',
    'KC242794': 'green',
    'KC242796': 'blue',
    'KC242799': 'blue',
    'KC242800': 'darkred',
    'KC242801': 'purple',
    'KC545389': 'lime',
    'KC545390': 'lime',
    'KC545391': 'lime',
    'KC545392': 'lime',
    'KC545393': 'silver',
    'KC545394': 'silver',
    'KC545395': 'silver',
    'KC545396': 'silver',
    'KC589025': 'lime',
    'KJ660346': 'red',
    'KJ660347': 'red',
    'KJ660348': 'red',
    'KM034555': 'red',
    'KM034557': 'red',
    'KM034560': 'red',
    'KM034562': 'red',
    'KM233039': 'red',
    'KM233050': 'red',
    'KM233053': 'red',
    'KM233057': 'red',
    'KM233063': 'red',
    'KM233070': 'red',
    'KM233072': 'red',
    'KM233096': 'red',
    'KM233097': 'red',
    'KM233099': 'red',
    'KM233103': 'red',
    'KM233109': 'red',
    'KM233110': 'red',
    'NC002549': 'purple',
}

for file_name in files:
    with open(path.join(BASE_PATH, file_name)) as f:
        data = f.readline()

    genome_name = file_name[:-11].replace('_', '')

    genome_color = colors[genome_name]

    meta.append(dict(file_name=file_name, genome_name=genome_name, color=genome_color))

with open(path.join(BASE_PATH, 'meta.json'), 'w') as f:
    json.dump(meta, f, indent=4)
