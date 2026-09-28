.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

GPU_BACKEND ?= cuda
SPFP ?= 0

# Compiler and architecture
ifeq ($(GPU_BACKEND),cuda)
GPU_CXX ?= nvcc
GPU_ARCH ?= sm_89
ARCH_FLAGS := -arch=$(GPU_ARCH)
GPU_LANG_FLAGS :=
BACKEND_DEFS := -DQGPU_BACKEND_CUDA

else ifeq ($(GPU_BACKEND),hip)
GPU_CXX ?= hipcc
GPU_ARCH ?=
ARCH_FLAGS := --offload-arch=$(GPU_ARCH)
GPU_LANG_FLAGS := -x hip
BACKEND_DEFS := -DQGPU_BACKEND_HIP -DHIP_ENABLE_WARP_SYNC_BUILTINS
export HIP_PLATFORM := amd

# HPC login nodes may have no GPU: require an explicit target.
ifneq ($(MAKECMDGOALS),clean)
ifeq ($(strip $(GPU_ARCH)),)
$(error HIP build requires GPU_ARCH, e.g. GPU_ARCH=gfx90a)
endif
endif

else
$(error Unsupported GPU_BACKEND=$(GPU_BACKEND); use cuda or hip)
endif

# Precision
ifeq ($(SPFP),1)
PRECISION := spfp
PRECISION_DEFS := -DQDYN_SPFP
else ifeq ($(SPFP),0)
PRECISION := double
PRECISION_DEFS :=
else
$(error SPFP must be 0 or 1)
endif

INCLUDES := -I./cuda/include -I./common/include -I./cpu/include -I.
BUILD_FLAGS := -O3 -std=c++17
DEFINES := $(BACKEND_DEFS) $(PRECISION_DEFS)

# Isolate backend, architecture and precision.
BUILD_DIR := .build/qgpu/$(GPU_BACKEND)/$(GPU_ARCH)/$(PRECISION)
TARGET := $(BUILD_DIR)/qdyn
BIN_DIR := ../../bin

CU_SRCS := main.cu $(wildcard cuda/src/*.cu)
CPP_SRCS := $(wildcard common/*.cpp common/src/*.cpp \
                       cpu/*.cpp cpu/src/*.cpp)

OBJS := $(addprefix $(BUILD_DIR)/,$(CU_SRCS:.cu=.o)) \
        $(addprefix $(BUILD_DIR)/,$(CPP_SRCS:.cpp=.o))
DEPS := $(OBJS:.o=.d)

DEP_FLAGS = -MMD -MF $(@:.o=.d) -MT $@

all: move

qdyn: $(TARGET)

$(TARGET): $(OBJS)
	$(GPU_CXX) $(BUILD_FLAGS) $(ARCH_FLAGS) $(LDFLAGS) \
		-o $@ $(OBJS) $(LDLIBS)

$(BUILD_DIR)/%.o: %.cu Makefile
	@mkdir -p $(@D)
	$(GPU_CXX) $(GPU_LANG_FLAGS) $(BUILD_FLAGS) $(ARCH_FLAGS) \
		$(INCLUDES) $(DEFINES) $(CPPFLAGS) $(CXXFLAGS) \
		$(DEP_FLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.cpp Makefile
	@mkdir -p $(@D)
	$(GPU_CXX) $(GPU_LANG_FLAGS) $(BUILD_FLAGS) $(ARCH_FLAGS) \
		$(INCLUDES) $(DEFINES) $(CPPFLAGS) $(CXXFLAGS) \
		$(DEP_FLAGS) -c $< -o $@

# Publish only after linking has completed.
# Keep the original executable inside its build directory.
move: qdyn
	@mkdir -p $(BIN_DIR)
	cp $(TARGET) $(BIN_DIR)/qdyn

double:
	$(MAKE) SPFP=0 all

spfp:
	$(MAKE) SPFP=1 all

clean:
	rm -rf .build/qgpu
	rm -f $(BIN_DIR)/qdyn

-include $(DEPS)

.PHONY: all qdyn move double spfp clean