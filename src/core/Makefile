CC = nvcc
SPFPFLAGS =
ifeq ($(QDYN_SPFP),1)
SPFPFLAGS += -DQDYN_SPFP
endif
CFLAGS = -O3 -std=c++17 -arch=sm_86 $(SPFPFLAGS) -I./cuda/include -I./common/include -I./cpu/include -I.
DEPFLAGS = -MMD -MF $(@:.o=.d)
BUILD_MODE = $(if $(filter 1,$(QDYN_SPFP)),spfp,dpfp)
OBJDIR = .build/$(BUILD_MODE)
TARGET = $(OBJDIR)/qdyn

# collect all .cu files except main.cu
SRCS = $(filter-out main.cu, $(wildcard *.cu cuda/src/*.cu))
CPPSRCS = $(wildcard common/*.cpp common/src/*.cpp cpu/*.cpp cpu/src/*.cpp)
MAIN_OBJ = $(OBJDIR)/main.o
OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.cu=.o)) $(addprefix $(OBJDIR)/,$(CPPSRCS:.cpp=.o))
DEPS = $(MAIN_OBJ:.o=.d) $(OBJS:.o=.d)

all: qdyn move

qdyn: $(TARGET)
	cp $< $@

$(TARGET): $(MAIN_OBJ) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJDIR)/%.o: %.cu
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.cpp
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

clean:
	rm -rf .build qdyn ../../bin/qdyn

move: $(TARGET)
	mkdir -p ../../bin
	cp $< ../../bin/qdyn

.PHONY: all qdyn clean move

-include $(DEPS)
