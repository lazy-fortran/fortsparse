#include "fortsparse_native.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#ifdef FS_TEST_COUNT_AMD
extern "C" int fsparse_test_general_amd_calls(void);
#endif
#define CHECK(x) do { if (!(x)) { std::fprintf(stderr,"line %d: %s\n",__LINE__,#x); std::exit(1); } } while(0)
int main() {
    using Z=std::complex<double>;
    static_assert(sizeof(Z)==2*sizeof(double),"C++ complex array representation");
    void *h=nullptr;
    int32_t cp[]={0,2,4},ri[]={0,1,0,1};
    Z v[]={{3,2},{1,-1},{2,-3},{5,1}};
    Z known[]={{1,2},{-2,1}},b[2],x[2];
    CHECK(fs_native_create(1,&h)==0);
    struct fs_native_memory_stats memory{};
    CHECK(fs_native_set_memory(h,0,16,FS_NATIVE_MEMORY_AUTO)==0);
    CHECK(fs_native_get_memory_stats(h,&memory)==0);
    CHECK(memory.host_limit_bytes==0 && memory.device_limit_bytes==16);
    CHECK(memory.host_current_bytes==0 && memory.device_current_bytes==0);
    CHECK(fs_native_set_refine(h,1)==0);
    CHECK(fs_native_factor_z(h,2,4,cp,ri,v)==0);
    char msg[160];
    int32_t nd=fs_native_ordering_available(FS_NATIVE_ORDER_ND);
#ifdef FS_TEST_EXPECT_METIS
    CHECK(nd==0);
#else
    CHECK(nd==3);
#endif
    if(nd==3) {
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_ND)==3);
        CHECK(fs_native_last_error(h,msg,sizeof msg)==3);
    }
    CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_DEFAULT)==0);
    CHECK(fs_native_last_error(h,msg,sizeof msg)==0 && msg[0]=='\0');
    int32_t amd=fs_native_ordering_available(FS_NATIVE_ORDER_AMD);
#ifdef FS_TEST_EXPECT_AMD
    CHECK(amd==0);
#else
    CHECK(amd==3);
#endif
    for(int policy=0;policy<3;++policy) {
        if(policy==1) {
#ifdef FS_TEST_COUNT_AMD
            int calls_before=fsparse_test_general_amd_calls();
#endif
            CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_AMD)==amd);
            if(amd==0) {
                CHECK(fs_native_solve_z(h,0,2,1,b,x)==4);
                CHECK(fs_native_factor_z(h,2,4,cp,ri,v)==0);
#ifdef FS_TEST_COUNT_AMD
                CHECK(fsparse_test_general_amd_calls()>calls_before);
#endif
                CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_AMD)==0);
            }
        } else if(policy==2) {
            CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_DEFAULT)==0);
            if(amd==0) CHECK(fs_native_solve_z(h,0,2,1,b,x)==4);
            CHECK(fs_native_refactor_z(h,2,4,cp,ri,v)==0);
        }
        for (int op=0;op<3;++op) {
            for(int i=0;i<2;++i) {
                b[i]=0;
                for(int j=0;j<2;++j) {
                    Z a=op ? v[i*2+j] : v[j*2+i];
                    b[i]+=(op==2 ? std::conj(a) : a)*known[j];
                }
            }
            CHECK(fs_native_solve_z(h,op,2,1,b,x)==0);
            CHECK(std::abs(x[0]-known[0])<1e-12 && std::abs(x[1]-known[1])<1e-12);
        }
    }
    CHECK(fs_native_free(&h)==0 && h==nullptr);
}
