import glob

files = sorted(glob.glob('results/md*/dualtop/output/energies.csv'))
files = files[::-1]
for d in files:
    dd = d.split('/')
    print('{} {}.json'.format(d,dd[1]))
