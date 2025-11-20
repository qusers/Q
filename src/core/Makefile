CC = nvcc
CFLAGS = -O3 -arch=sm_89 -I./cuda/include -I.

# collect all .cu files except main.cu
SRCS = $(filter-out main.cu, $(wildcard *.cu cuda/src/*.cu))
OBJS = $(SRCS:.cu=.o)

all: qdyn move

qdyn: main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cu
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o cuda/src/*.o ../../bin/qdyn

move:
	mkdir -p ../../bin
	mv qdyn ../../bin/
