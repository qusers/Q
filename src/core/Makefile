CC = nvcc
CFLAGS = -O3 -std=c++17 -arch=sm_86 -I./cuda/include -I./common/include -I./cpu/include -I.
DEPFLAGS = -MMD -MF $(@:.o=.d)
QDYN_SPFP ?= 0

ifneq ($(filter 1 true TRUE yes YES on ON,$(QDYN_SPFP)),)
CFLAGS += -DQDYN_SPFP
PRECISION = spfp
else
PRECISION = dpfp
endif

BUILD_DIR = build/$(PRECISION)

# collect all .cu files except main.cu
SRCS = $(filter-out main.cu, $(wildcard *.cu cuda/src/*.cu))
CPPSRCS = $(wildcard common/*.cpp common/src/*.cpp cpu/*.cpp cpu/src/*.cpp)
MAIN_OBJ = $(BUILD_DIR)/main.cu.o
OBJS = $(addprefix $(BUILD_DIR)/,$(SRCS:.cu=.cu.o)) $(addprefix $(BUILD_DIR)/,$(CPPSRCS:.cpp=.cpp.o))
DEPS = $(MAIN_OBJ:.o=.d) $(OBJS:.o=.d)

all: qdyn move

qdyn: $(MAIN_OBJ) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/%.cu.o: %.cu
	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

clean:
	rm -rf build
	rm -f *.o *.d cuda/src/*.o cuda/src/*.d common/*.o common/*.d common/src/*.o common/src/*.d cpu/*.o cpu/*.d cpu/src/*.o cpu/src/*.d ../../bin/qdyn

move:
	mkdir -p ../../bin
	mv qdyn ../../bin/

-include $(DEPS)
