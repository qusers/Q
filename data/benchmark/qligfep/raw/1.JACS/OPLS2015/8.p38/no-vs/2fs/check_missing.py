l1 = []
l2 = []

with open('protein.txt', 'r') as infile:
	for line in infile:			
		line = line.split()
		l1.append(line[0])
infile.close()

for f in ['water.txt']:
	with open(f, 'r') as infile:
		for line in infile:
			line = line.split()
			l2.append(line[0])
infile.close()

list1 = [item for item in l2 if item not in l1] #in water but not in protein
list2 = [item for item in l1 if item not in l2] #in protein but not in water
print(list2)

with open('protein.txt', 'r') as infile:
	with open('super_protein.txt', 'w') as outfile1:
		for line in infile:
			line = line.split()
			if line[0] not in list2:
				outfile1.write('{}\t{}\t{}\t{}\n'.format(line[0], line[1], line[2], line[3]))
			else:
				continue
infile.close()
outfile1.close()

with open('water.txt', 'r') as infile:
        with open('super_water.txt', 'w') as outfile2:
                for line in infile:
                        line = line.split()
                        if line[0] not in list1:
                                outfile2.write('{}\t{}\t{}\t{}\n'.format(line[0], line[1], line[2], line[3]))
                        else:
                                continue
infile.close()
outfile2.close()


