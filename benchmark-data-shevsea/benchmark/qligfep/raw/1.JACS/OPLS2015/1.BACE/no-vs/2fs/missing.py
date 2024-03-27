list1 = []
with open('protein.txt', 'r') as infile:
    for line in infile:
        line = line.split()
        list1.append(line[0])

list2 = []
with open('water.txt', 'r') as infile2:
    for line in infile2:
        line = line.split()
        list2.append(line[0])

nl = [i for i in list1 if i not in list2]

for i in nl:
    print('{}\n'.format(i))

print(set([x for x in list2 if list2.count(x) > 1]))
