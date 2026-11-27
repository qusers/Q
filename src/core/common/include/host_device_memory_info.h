#pragma once
class HostDeviceMemoryInfo {
   public:
    long long total_cpu_memory;
    long long total_gpu_memory;

    static HostDeviceMemoryInfo& Instance() {
        static HostDeviceMemoryInfo instance;
        return instance;
    }

   protected:
    HostDeviceMemoryInfo() {
        total_cpu_memory = total_gpu_memory = 0;
    }
};