import os
from os import path
import re
import json


BASE_PATH = './Genomes/corona'
EXT = '.fasta'


files = [f for f in os.listdir(BASE_PATH) if path.isfile(path.join(BASE_PATH, f)) and f.endswith(EXT)]
meta = []
colors = {
    'Murine hepatitis virus strain 2': 'green',
    'Murine hepatitis virus strain Penn 97-1': 'green',
    'Murine hepatitis virus strain ML-10': 'green',
    'Bovine coronavirus strain Quebec': 'green',
    'Human coronavirus 229E': 'red',
    'Porcine epidemic diarrhea virus strain CV777': 'red',
    'Bovine coronavirus isolate BCoV-LUN': 'green',
    'SARS coronavirus BJ01': 'blue',
    'SARS coronavirus HKU-39849': 'blue',
    'SARS coronavirus CUHK-W1': 'blue',
    'SARS coronavirus Urbani': 'blue',
    'SARS coronavirus CUHK-Su10': 'blue',
    'SARS coronavirus Sin2500': 'blue',
    'SARS coronavirus Sin2677': 'blue',
    'SARS coronavirus Sin2679': 'blue',
    'SARS coronavirus Sin2748': 'blue',
    'SARS coronavirus Sin2774': 'blue',
    'SARS coronavirus TW1': 'blue',
    'SARS coronavirus ZJ01': 'blue',
    'Human coronavirus OC43': 'green',
    'SARS coronavirus civet007': 'blue',
    'SARS coronavirus civet010': 'blue',
    'Turkey coronavirus isolate MG10': 'olive',
    'Avian infectious bronchitis virus': 'olive',
    'O\'nyong-nyong virus': 'purple',
    'Ross River virus': 'purple',
    'Cell fusing agent virus strain Galveston': 'purple',
    'Mouse hepatitis virus strain MHV-A59 C12 mutant': 'green',
    'Bovine coronavirus': 'green',
    'Hepatitis C virus genotype 1': 'purple',
    'SARS coronavirus': 'blue',
    'Human Coronavirus NL63': 'red',
    'Human coronavirus HKU1': 'darkred',
    'Bovine coronavirus strain Mebus': 'green',
}
small_name = {
    'Murine hepatitis virus strain 2': '2 MHV2',
    'Murine hepatitis virus strain Penn 97-1': '2 MHVP',
    'Murine hepatitis virus strain ML-10': '2 MHVM',
    'Bovine coronavirus strain Quebec': '2 BCoVQ',
    'Human coronavirus 229E': '1 HCoV-229E',
    'Porcine epidemic diarrhea virus strain CV777': '1 PEDV',
    'Bovine coronavirus isolate BCoV-LUN': '2 BCoVL',
    'SARS coronavirus BJ01': '4 BJ01',
    'SARS coronavirus HKU-39849': '4 HKU-39849',
    'SARS coronavirus CUHK-W1': '4 CUHK-W1',
    'SARS coronavirus Urbani': '4 Urbani',
    'SARS coronavirus CUHK-Su10': '4 CUHK-Su10',
    'SARS coronavirus Sin2500': '4 Sin2500',
    'SARS coronavirus Sin2677': '4 Sin2677',
    'SARS coronavirus Sin2679': '4 Sin2679',
    'SARS coronavirus Sin2748': '4 Sin2748',
    'SARS coronavirus Sin2774': '4 Sin2774',
    'SARS coronavirus TW1': '4 TW1',
    'SARS coronavirus ZJ01': '4 ZJ01',
    'Human coronavirus OC43': '2 HCoV-OC43',
    'SARS coronavirus civet007': '4 civet007',
    'SARS coronavirus civet010': '4 civet010',
    'Turkey coronavirus isolate MG10': '3 TCoV',
    'Avian infectious bronchitis virus': '3 IBV',
    'O\'nyong-nyong virus': 'out NyongT',
    'Ross River virus': 'out RossT',
    'Cell fusing agent virus strain Galveston': 'out CellF',
    'Mouse hepatitis virus strain MHV-A59 C12 mutant': '2 MHV',
    'Bovine coronavirus': '2 BCoV',
    'Hepatitis C virus genotype 1': 'out HepaCF',
    'SARS coronavirus': '4 TOR2',
    'Human Coronavirus NL63': '1 HCoV-NL63',
    'Human coronavirus HKU1': '5 HCoV-HKU1',
    'Bovine coronavirus strain Mebus': '2 BCoVM',
}

for file_name in files:
    with open(path.join(BASE_PATH, file_name)) as f:
        data = f.readline()

    genome_name = re.search(r'(?<=\.[0-9]\s).+(?=,\s)', data)[0]

    genome_color = colors[genome_name]

    meta.append(dict(file_name=file_name, genome_name=small_name[genome_name], color=genome_color))

with open(path.join(BASE_PATH, 'meta.json'), 'w') as f:
    json.dump(meta, f, indent=4)
