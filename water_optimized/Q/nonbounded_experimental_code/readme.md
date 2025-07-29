# Non-bounded kernel code
These codes use different features to calculate non-bounded forces.
1. `old.cu`: It's the old version non-bounded kernel.
2. `v1.cu`: It merges three kernels into one kernel and makes one thread execute more tile atoms, without using any features, like shared memory.
3. `v2.cu`: It uses some features like shared memory, warp shuffle, and instruction optimization to optimize the `v1.cu`.
4. `v3.cu`: Use symmetry-aware folding algorithm to change the shape of the kernels. In `v3`, there are two parameters to control the size of atom pairs for one thread to calculate. `Block_x` and `Block_y`. It means it can not only calculate a column of values. It can also control the number of rows.
5. `v4.cu`: In this version, it will make `Block_x = 1` to decrease the shared memory usage. It can increase the occupancy of SM. But it seems doesn't have some optimization. 