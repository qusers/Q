import json

block = 0
Q_energies = {}
velocities = []

with open('eq1.log') as infile:
    for line in infile:
        if len(line.strip()) < 2:
            continue        
        
        if 'k (vk)' in line:
            line = line.split()
            velocities.append(line[3])
            
        if 'Energy summary at step      0' in line:
            block = 1
            continue
        if 'FINAL' in line:
            block = 0
            break
            
        if block == 1:
            line = line.split()
            if line[0] == 'solute':
                Q_energies['solute'] = line[1:]
                
            if line[0] == 'solvent':
                Q_energies['solvent'] = line[1:]
                                
            if line[0] == 'solute-solvent':
                Q_energies['solute-solvent'] = line[1:]        
                                
            if line[0] == 'restraints':
                Q_energies['restraints'] = line[1:]
                                                
            if line[0] == 'SUM':
                Q_energies['SUM'] = line[1:]

data = json.dumps(Q_energies, indent=1)
with open('Q_data.json', 'w') as outfile:
    outfile.write(data)

with open('velocities.csv','w') as outfile:
    outfile.write('{}\n'.format(int(len(velocities)/3)))
    outfile.write('0\n')
    for i, v in enumerate(velocities):
        if (i + 1) % 3 != 0:
            outfile.write('{};'.format(v))
        
        else:
            outfile.write('{}\n'.format(v))
