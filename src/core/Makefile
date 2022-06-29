CC = nvcc
CFLAGS = -g

all: qdyn move

qdyn: system.o main.o utils.o parse.o restraints.o bonded.o nonbonded.o solvent.o qatoms.o shake.o
	$(CC) $(CFLAGS) -o qdyn main.o system.o utils.o parse.o restraints.o bonded.o nonbonded.o solvent.o qatoms.o shake.o

main.o: main.cu system.h
	$(CC) $(CFLAGS) -c main.cu

system.o: system.cu system.h utils.h parse.h restraints.h
	$(CC) $(CFLAGS) -c system.cu

utils.o: utils.cu utils.h
	$(CC) $(CFLAGS) -c utils.cu

parse.o: parse.cu parse.h system.h
	$(CC) $(CFLAGS) -c parse.cu

restraints.o: restraints.cu restraints.h system.h
	$(CC) $(CFLAGS) -c restraints.cu

bonded.o: bonded.cu bonded.h system.h
	$(CC) $(CFLAGS) -c bonded.cu

nonbonded.o: nonbonded.cu nonbonded.h system.h
	$(CC) $(CFLAGS) -c nonbonded.cu

solvent.o: solvent.cu solvent.h system.h
	$(CC) $(CFLAGS) -c solvent.cu

qatoms.o: qatoms.cu qatoms.h system.h
	$(CC) $(CFLAGS) -c qatoms.cu

shake.o: shake.cu system.h
	$(CC) $(CFLAGS) -c shake.cu

clean:
	rm -f *.o ../../bin/qdyn

move:
	mkdir -p ../../bin
	mv qdyn ../../bin/
