CC = nvcc
CFLAGS = -g

all: qdyn

qdyn: system.o main.o utils.o parse.o restraints.o bonded.o nonbonded.o solvent.o
	$(CC) $(CFLAGS) -o qdyn main.o system.o utils.o parse.o restraints.o bonded.o nonbonded.o solvent.o

main.o: main.cu system.h
	$(CC) $(CFLAGS) -c main.cu

system.o: system.cu system.h utils.h parse.h restraints.h
	$(CC) $(CFLAGS) -c system.cu

utils.o: utils.cu utils.h
	$(CC) $(CFLAGS) -c utils.cu

parse.o: parse.cu parse.h
	$(CC) $(CFLAGS) -c parse.cu

restraints.o: restraints.cu restraints.h system.h
	$(CC) $(CFLAGS) -c restraints.cu

bonded.o: bonded.cu bonded.h system.h
	$(CC) $(CFLAGS) -c bonded.cu

nonbonded.o: nonbonded.cu nonbonded.h system.h
	$(CC) $(CFLAGS) -c nonbonded.cu

solvent.o: solvent.cu solvent.h system.h
	$(CC) $(CFLAGS) -c solvent.cu

clean:
	rm -f *.o qdyn
