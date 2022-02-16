import glob
import os
import numpy as np

results = {'old':{}, 'new':{}}
target_edges = {}
for i, target in enumerate(['cdk8', 'cmet', 'eg5', 'hif2a', 'pfkfb3', 'shp2', 'syk', 'tnks2'], start=1):
    target_edges[target] = []
    with open('maps/{}_MFP.txt'.format(target), 'r') as infile:
        for line in infile:
            line = line.split()
            target_edges[target].append((line[0], line[1]))

    results['old'][target] = {}
    for leg in ['protein', 'water']:
        results['old'][target][leg] = {}
        with open('{}/no-vs/2fs/{}.txt'.format(target, leg), 'r') as infile:
            for line in infile:
                line = line.split()
                if line[0] == 'FEP':
                    continue
                else:
                    try:
                        fep = line[0].split('FEP_')[1]
                        l1, l2 = fep.split('_')[0], fep.split('_')[1]
                    except:
                        fep = line[0].split('FEP_')[1]
                        l1, l2 = fep.split('-')[0], fep.split('-')[1]
                    if (l1, l2) in target_edges[target]:
                        fep = line[0]
                        if target == 'CDK2':
                            fep = line[0].replace('-', '_')
                        results['old'][target][leg][fep] = [float(line[1]), float(line[2]), int(line[3])]
                    elif (l2, l1) in target_edges[target]:
                        results['old'][target][leg]['FEP_{}_{}'.format(l2, l1)] = [-1*float(line[1]), float(line[2]), int(line[3])]
    results['new'][target] = {}
    for leg in ['protein', 'water']:
        results['new'][target][leg] = {}
        with open('{}/benchmark/2fs/{}.txt'.format(target, leg), 'r') as infile:
            for line in infile:
                line = line.split()
                if line[0] == 'FEP':
                    continue
                else:
                    try:
                        fep = line[0].split('FEP_')[1]
                        l1, l2 = fep.split('_')[0], fep.split('_')[1]
                    except:
                        fep = line[0].split('FEP_')[1]
                        l1, l2 = fep.split('-')[0], fep.split('-')[1]
                    if (l1, l2) in target_edges[target]:
                        results['new'][target][leg][line[0]] = [float(line[1]), float(line[2]), int(line[4])]
                    elif (l2, l1) in target_edges[target]: 
                        results['new'][target][leg]['FEP_{}_{}'.format(l2, l1)] = [-1*float(line[1]), float(line[2]), int(line[4])]

results['merged'] = {}
for i, target in enumerate(['cdk8', 'cmet', 'eg5', 'hif2a', 'pfkfb3', 'shp2', 'syk', 'tnks2'], start=1):
    results['merged'][target] = {}
    for leg in ['protein', 'water']:
        results['merged'][target][leg] = {}
        n, o = results['new'][target][leg], results['old'][target][leg]
        for edge in target_edges[target]:
            try:
                try:
                    fep = 'FEP_{}_{}'.format(edge[0], edge[1])
                    ddG = sum([n[fep][0], o[fep][0]])/2
                    sem = np.sqrt(sum([n[fep][1]**2, o[fep][1]**2]))
                    crashes = sum([n[fep][2], o[fep][2]])
                    std = sem*np.sqrt(20 - crashes)
                    results['merged'][target][leg][fep] = [ddG, sem, std, crashes]
                except:
                    fep = 'FEP_{}_{}'.format(edge[1], edge[0])
                    ddG = sum([n[fep][0], o[fep][0]])/2
                    sem = np.sqrt(sum([n[fep][1]**2, o[fep][1]**2]))
                    crashes = sum([n[fep][2], o[fep][2]])
                    std = sem*np.sqrt(20 - crashes)
                    results['merged'][target][leg][fep] = [ddG, sem, std, crashes]
            except:
                continue

        merged_dir = os.path.join(os.getcwd(), target, 'merged', '2fs')       
        if not os.path.exists(merged_dir):
            os.mkdir(os.path.join(os.getcwd(), target, 'merged'))
            os.mkdir(os.path.join(os.getcwd(), target, 'merged', '2fs'))
        
        with open('{}/{}.txt'.format(merged_dir, leg), 'w') as outfile:
            outfile.write('{:<70}{:<10}{:<10}{:<10}{:<10}\n'.format('FEP','dG','sem', 'std', 'crashes'))
            for fep, v in results['merged'][target][leg].items():
                outfile.write('{:<70}{:<10.3}{:<10.3}{:<10.3}{:<10}\n'.format(fep, v[0], v[1], v[2], v[3]))
