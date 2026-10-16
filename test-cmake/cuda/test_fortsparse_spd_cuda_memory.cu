#include "fsparse_spd_cuda.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>

namespace {
std::map<void*,size_t> allocations;
int calls, fail_call;
size_t live, high;
void require(bool ok) { if (!ok) { std::fputs("CUDA payload oracle failed\n",stderr); std::exit(1); } }
void check(int code) { if (code) { std::fprintf(stderr,"CUDA payload status %d\n",code); std::exit(1); } }
void reset(int failure=0) { require(allocations.empty()); calls=0;fail_call=failure;live=high=0; }
fsparse_spd_cuda_memory stats(void *ctx) {
    fsparse_spd_cuda_memory m={};check(fsparse_spd_cuda_memory_stats(ctx,&m));
    require(size_t(m.current_bytes)==live && size_t(m.peak_bytes)==high);
    require(fsparse_spd_cuda_allocated_bytes(ctx)==live);return m;
}
const fsparse_spd_cuda_level_panel panels[]={{0,1,2,2,0},{2,1,1,1,1}};
const int maps[]={1},tail[]={2};
const fsparse_spd_cuda_schur_task task={1,1,2,0,0,2,2,1,1,1,1,1};
const fsparse_spd_cuda_level_task level_task={task,0,1};
const size_t ptr[]={0,1};
int plan(void *ctx,int kind) {
    if(kind==0) return fsparse_spd_cuda_set_schur_batches(ctx,&task,1,ptr,1);
    if(kind==1) return fsparse_spd_cuda_set_levels(ctx,panels,2,&level_task,1,2);
    return fsparse_spd_cuda_set_solve_plan(ctx,2,panels,2,tail,1,2);
}
void numerical(void *ctx) {
    // Original A=[[4,2],[2,10]], x=[1,2], b=[8,22]. Independent known-x oracle.
    const double a[]={4,2,10}, b[]={8,22}; double x[2]={};
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_upload(ctx,0,a,3));
    for(size_t i=0;i<2;++i) {
        check(fsparse_spd_cuda_factor_level(ctx,i));
        check(fsparse_spd_cuda_update_level(ctx,i));
    }
    int info=-1;check(fsparse_spd_cuda_factor_info(ctx,&info));require(info==0);
    check(fsparse_spd_cuda_upload(ctx,8,b,2));check(fsparse_spd_cuda_solve(ctx,8,2,1));
    check(fsparse_spd_cuda_download(ctx,8,x,2));check(fsparse_spd_cuda_synchronize(ctx));
    require(std::abs(x[0]-1)<1e-13 && std::abs(x[1]-2)<1e-13);
    require(std::abs(4*x[0]+2*x[1]-8)<1e-12 && std::abs(2*x[0]+10*x[1]-22)<1e-12);
}
int workflow(size_t cap,int failure=0) {
    reset(failure);void *ctx=nullptr;fsparse_spd_cuda_memory report={};
    int code=fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,cap,&report);
    if(!code) stats(ctx);
    if(!code) { code=fsparse_spd_cuda_set_maps(ctx,maps,1);report=stats(ctx); }
    for(int pass=0;pass<2 && !code;++pass)
        for(int k=0;k<3 && !code;++k) { code=plan(ctx,k);report=stats(ctx); }
    // Exercise the explicit workspace last: fault injection then avoids paying
    // repeated opaque cuBLAS handle initialization for unrelated plan faults.
    if(!code) { code=fsparse_spd_cuda_prepare_blas(ctx);report=stats(ctx); }
    if(!code) numerical(ctx);
    if(code==-2) require(report.required_bytes>report.available_bytes);
    if(ctx) check(fsparse_spd_cuda_destroy(ctx));
    else require(report.current_bytes==0 && size_t(report.peak_bytes)==high);
    require(allocations.empty() && live==0);
    return code;
}
void replacement(int kind) {
    // Learn request sizes independently via wrapped CUDA allocations, then use
    // the exact old+new peak minus one. The old plan must survive refusal.
    reset();void *ctx=nullptr;fsparse_spd_cuda_memory report={};
    check(fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,0,&report));
    check(fsparse_spd_cuda_set_maps(ctx,maps,1));check(plan(ctx,kind));
    check(plan(ctx,kind));const size_t peak=high;
    check(fsparse_spd_cuda_destroy(ctx));require(allocations.empty());
    reset();check(fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,peak-1,&report));
    check(fsparse_spd_cuda_set_maps(ctx,maps,1));check(plan(ctx,kind));
    const auto old=allocations;const size_t before=live;
    require(plan(ctx,kind)==-2);report=stats(ctx);
    require(allocations==old && live==before && report.required_bytes>report.available_bytes);
    // Clear and rebuild under the same cap proves recovery without increasing it.
    if(kind==0) check(fsparse_spd_cuda_set_schur_batches(ctx,nullptr,0,ptr,0));
    else if(kind==1) check(fsparse_spd_cuda_set_levels(ctx,nullptr,0,nullptr,0,0));
    else { check(fsparse_spd_cuda_set_maps(ctx,maps,1)); }
    check(plan(ctx,kind));stats(ctx);
    check(fsparse_spd_cuda_destroy(ctx));require(allocations.empty() && live==0);
}
}
extern "C" cudaError_t __real_cudaMalloc(void**,size_t);
extern "C" cudaError_t __real_cudaFree(void*);
extern "C" cudaError_t __wrap_cudaMalloc(void **p,size_t bytes) {
    ++calls;if(calls==fail_call) { *p=nullptr;return cudaErrorMemoryAllocation; }
    cudaError_t code=__real_cudaMalloc(p,bytes);
    if(code==cudaSuccess) { allocations[*p]=bytes;live+=bytes;high=std::max(high,live); }
    return code;
}
extern "C" cudaError_t __wrap_cudaFree(void *p) {
    cudaError_t code=__real_cudaFree(p);
    if(code==cudaSuccess && p) { auto i=allocations.find(p);require(i!=allocations.end());live-=i->second;allocations.erase(i); }
    return code;
}
int main() {
    check(workflow(0));const size_t exact=high;const int count=calls;
    check(workflow(exact));require(high==exact);
    require(workflow(exact-1)==-2);
    for(int i=1;i<=count;++i) {
        require(workflow(0,i)==int(cudaErrorMemoryAllocation));require(calls==i);
    }
    for(int k=0;k<3;++k) replacement(k);
    // Initial arena and info refusals return no context or leaked allocation.
    reset();void *ctx=nullptr;fsparse_spd_cuda_memory report={};
    require(fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,1023,&report)==-2);
    require(!ctx && report.current_bytes==0 && report.peak_bytes==0 && live==0);
    require(fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,1027,&report)==-2);
    require(!ctx && report.current_bytes==0 && report.peak_bytes==1024 && live==0);
    require(fsparse_spd_cuda_create_deferred_limited(&ctx,0,std::numeric_limits<size_t>::max(),0,&report)==-1);
    reset();check(fsparse_spd_cuda_create_deferred_limited(&ctx,0,128,0,&report));
    const int before=calls;
    const size_t too_many=size_t(INT64_MAX)/sizeof(int)+1;
    if(too_many<=std::numeric_limits<size_t>::max()/sizeof(int)) {
        require(fsparse_spd_cuda_set_maps(ctx,maps,too_many)==-3);
        require(calls==before);stats(ctx);
    }
    check(fsparse_spd_cuda_destroy(ctx));require(live==0 && allocations.empty());
    std::printf("CUDA exact payload caps, replacement overlap, %d allocation faults and original-A oracle pass\n",count);
}
