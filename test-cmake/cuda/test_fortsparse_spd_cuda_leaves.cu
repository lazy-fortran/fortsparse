#include "fortsparse_spd_leaves.cuh"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

__global__ void evaluate(double* values)
{
    const double a = 4.0, b = -2.0, c = 0.75;
    const double da = -0.3, db = 0.2, dc = -0.1, seed = 1.25;
    spd_root_jvp(a, da, values, values+1);
    spd_root_vjp(a, seed, values+2, values+3);
    spd_divide_jvp(a, b, da, db, values+4, values+5);
    spd_divide_vjp(a, b, seed, values+6, values+7, values+8);
    spd_update_jvp(a, b, c, da, db, dc, values+9, values+10);
    spd_update_vjp(a, b, c, seed, values+11, values+12, values+13, values+14);
}

int main()
{
    double *device = nullptr;
    double host[15];
    const double a = 4.0, b = -2.0, c = 0.75;
    const double da = -0.3, db = 0.2, dc = -0.1, seed = 1.25;
    const double h = 1.0e-5;
    if (cudaMalloc(&device, sizeof(host)) != cudaSuccess) return 1;
    evaluate<<<1,1>>>(device);
    cudaError_t code = cudaGetLastError();
    if (code == cudaSuccess)
        code = cudaMemcpy(host, device, sizeof(host), cudaMemcpyDeviceToHost);
    cudaFree(device);
    if (code != cudaSuccess) {
        std::fprintf(stderr, "CUDA generated leaf test: %s\n", cudaGetErrorString(code));
        return 1;
    }
    const double expected[15] = {
        std::sqrt(a), (std::sqrt(a+h*da)-std::sqrt(a-h*da))/(2*h),
        std::sqrt(a), seed/(2*std::sqrt(a)),
        a/b, ((a+h*da)/(b+h*db)-(a-h*da)/(b-h*db))/(2*h),
        a/b, seed/b, -seed*a/(b*b),
        a-b*c, ((a+h*da)-(b+h*db)*(c+h*dc)
                  -((a-h*da)-(b-h*db)*(c-h*dc)))/(2*h),
        a-b*c, seed, -seed*c, -seed*b
    };
    for (int i = 0; i < 15; ++i) {
        if (!std::isfinite(host[i]) || std::abs(host[i]-expected[i]) > 2.0e-9) {
            std::fprintf(stderr, "CUDA generated leaf mismatch %d: %.17g vs %.17g\n",
                         i, host[i], expected[i]);
            return 1;
        }
    }
    if (std::abs(seed*host[1]-host[3]*da) > 2.0e-14 ||
        std::abs(seed*host[5]-host[7]*da-host[8]*db) > 2.0e-14 ||
        std::abs(seed*host[10]-host[12]*da-host[13]*db-host[14]*dc) > 2.0e-14)
        return 1;
    std::puts("CUDA SPD generated values, JVP/VJP and adjoint pairings pass");
    return 0;
}
