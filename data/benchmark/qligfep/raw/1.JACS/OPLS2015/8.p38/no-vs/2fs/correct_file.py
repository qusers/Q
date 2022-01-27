feps = []
with open('water.txt', 'r') as infile:
	for line in infile:
		line = line.split()
		if line[0] == 'FEP':
			continue
		else:
			feps.append(line[0])
infile.close()
for i in feps:
	if feps.count(i) > 1:
		print(i)
