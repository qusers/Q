import json

refdata = {}

#Q5toGPUmask = {
#                ''el       vdW      bond     angle   torsion  improper
#}

with open('Q5_data/Q_data.json') as infile:
    Q_data = json.load(infile)
    
with open('refdata.txt') as infile:
    for line in infile:
        line = line.strip()
        line = line.split('=')
        if len(line) > 1:
            refdata[line[0]] = line[1]
        
print(Q_data)
