import os
from os import path
import re
import json


BASE_PATH = './Genomes/Influenza A viruses'
EXT = '.fasta'


files = [f for f in os.listdir(BASE_PATH) if path.isfile(path.join(BASE_PATH, f)) and f.endswith(EXT)]
meta = []
colors = {
    'H5N1': 'red',
    'H1N1': 'blue',
    'H2N2': 'magenta',
    'H7N3': 'green',
    'H7N9': 'black',
}

for file_name in files:
    with open(path.join(BASE_PATH, file_name)) as f:
        data = f.readline()

    genome_name = re.search(r'\([^\(]*\(([^\)]*)\)\)', data)[0][1:-1]
    genome_name = genome_name.replace(' ', '')

    genome_color = 'white'
    for pattern, color in colors.items():
        if pattern in genome_name:
            genome_color = color
            break

    meta.append(dict(file_name=file_name, genome_name=genome_name, color=genome_color))

with open(path.join(BASE_PATH, 'meta.json'), 'w') as f:
    json.dump(meta, f, indent=4)
