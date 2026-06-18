#pragma once
#include "cpu_shake.h"

class CudaShakeV2 : public CpuShake {
    public:
    void apply(Context &ctx) override;
};