import os
import glob

cwd = os.getcwd()
missing = []
for leg in ['protein', 'water']:
    f = glob.glob('*{}.txt'.format(leg))[0]
    with open(f, 'r') as infile:
        for line in infile:
            frg = line.split()
            if frg[0] == 'FEP':
                continue
            if int(frg[3]) > 2:
                missing.append(frg[0][4:])

missing = list(set(missing))

for i in missing:
    print(i.replace('_', ' '))
