#include "fortsparse_native.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#define CHECK(x) do { if (!(x)) { std::fprintf(stderr,"line %d: %s\n",__LINE__,#x); std::exit(1); } } while(0)
int main() {
    using Z=std::complex<double>;
    static_assert(sizeof(Z)==2*sizeof(double),"C++ complex array representation");
    void *h=nullptr;
    int32_t cp[]={0,2,4},ri[]={0,1,0,1};
    Z v[]={{3,2},{1,-1},{2,-3},{5,1}};
    Z known[]={{1,2},{-2,1}},b[2],x[2];
    CHECK(fs_native_create(1,&h)==0);
    CHECK(fs_native_set_refine(h,1)==0);
    CHECK(fs_native_factor_z(h,2,4,cp,ri,v)==0);
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
    CHECK(fs_native_free(&h)==0 && h==nullptr);
}
