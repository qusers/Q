CC = nvcc
CFLAGS = -O3 -std=c++17 -arch=sm_89 -I./cuda/include -I./common/include -I./cpu/include -I.

# collect all .cu files except main.cu
SRCS = $(filter-out main.cu, $(wildcard *.cu cuda/src/*.cu))
CPPSRCS = $(wildcard common/*.cpp common/src/*.cpp cpu/*.cpp cpu/src/*.cpp)
OBJS = $(SRCS:.cu=.o) $(CPPSRCS:.cpp=.o)

all: qdyn move

qdyn: main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cu
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o cuda/src/*.o common/*.o common/src/*.o cpu/*.o cpu/src/*.o ../../bin/qdyn

move:
	mkdir -p ../../bin
	mv qdyn ../../bin/
