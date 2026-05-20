CC = nvcc
CFLAGS = -O3 -std=c++17 -arch=sm_86 -I./cuda/include -I./common/include -I./cpu/include -I.
DEPFLAGS = -MMD -MF $(@:.o=.d)
SPFP ?= 0

ifeq ($(SPFP),1)
CFLAGS += -DQDYN_SPFP
endif

# collect all .cu files except main.cu
SRCS = $(filter-out main.cu, $(wildcard *.cu cuda/src/*.cu))
CPPSRCS = $(wildcard common/*.cpp common/src/*.cpp cpu/*.cpp cpu/src/*.cpp)
MAIN_OBJ = main.o
OBJS = $(SRCS:.cu=.o) $(CPPSRCS:.cpp=.o)
DEPS = $(MAIN_OBJ:.o=.d) $(OBJS:.o=.d)

all: qdyn move

double:
	$(MAKE) clean
	$(MAKE) SPFP=0 all

spfp:
	$(MAKE) clean
	$(MAKE) SPFP=1 all

qdyn: $(MAIN_OBJ) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cu
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

clean:
	rm -f *.o *.d cuda/src/*.o cuda/src/*.d common/*.o common/*.d common/src/*.o common/src/*.d cpu/*.o cpu/*.d cpu/src/*.o cpu/src/*.d ../../bin/qdyn

move:
	mkdir -p ../../bin
	mv qdyn ../../bin/

-include $(DEPS)

.PHONY: all clean move double spfp
