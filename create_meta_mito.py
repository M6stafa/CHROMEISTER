import os
from os import path
import re
import json


BASE_PATH = './Genomes/mito'
EXT = '.fasta'


files = [f for f in os.listdir(BASE_PATH) if path.isfile(path.join(BASE_PATH, f)) and f.endswith(EXT)]
meta = []
colors = {
    'Sheep': 'magenta',
    'BrownBear': 'blue',
    'PolarBear': 'blue',
    'Goat': 'magenta',
    'Dormouse': 'black',
    'Rabbit': 'purple',
    'Pig': 'magenta',
    'Squirrel': 'black',
    'Buffalo': 'magenta',
    'VerMonkey': 'red',
    'ComChim': 'red',
    'Gorilla': 'red',
    'BorOran': 'red',
    'PigChim': 'red',
    'BlackBear': 'blue',
    'GiantPanda': 'blue',
    'Leopard': 'blue',
    'Tiger': 'blue',
    'Wolf': 'blue',
    'FinWhale': 'green',
    'horse': 'cyan',
    'donkey': 'cyan',
    'SumOran': 'red',
    'MacacaThibet': 'red',
    'BowheadWhale': 'green',
    'GrayWhale': 'green',
    'IndusRiverDolphin': 'green',
    'NorthPacificWhale': 'green',
    'Chiru': 'magenta',
    'CommonWarthog': 'magenta',
    'TaiwanSerow': 'magenta',
    'Cat': 'blue',
    'Dog': 'blue',
    'Cow': 'magenta',
    'Human': 'red',
    'BlueWhale': 'green',
    'Hedgehog': 'gray',
    'IndianRhino': 'cyan',
    'Gibbon': 'red',
    'WhiteRhino': 'cyan',
    'Baboon': 'red',
}

for file_name in files:
    with open(path.join(BASE_PATH, file_name)) as f:
        data = f.readline()

    genome_name = re.search(r'(?<=[0-9])[a-zA-Z]+(?=.fasta)', file_name)[0]
    genome_name = genome_name.replace(' ', '')

    genome_color = colors[genome_name]

    meta.append(dict(file_name=file_name, genome_name=genome_name, color=genome_color))

with open(path.join(BASE_PATH, 'meta.json'), 'w') as f:
    json.dump(meta, f, indent=4)
