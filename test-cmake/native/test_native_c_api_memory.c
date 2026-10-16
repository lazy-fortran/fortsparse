#include "fortsparse_native.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);exit(1);} } while(0)
int fsparse_native_amd_bounded_available(void);
void fsparse_spd_probe_begin(int64_t);
void fsparse_spd_probe_end(int64_t*,int64_t*,int64_t*);
static int nd_calls;
int __real_fsparse_native_order_general(int,int,const int*,const int*,int*,int*);
int __wrap_fsparse_native_order_general(int n,int nz,const int *ptr,const int *row,int *p,int *ip) {
    ++nd_calls;return __real_fsparse_native_order_general(n,nz,ptr,row,p,ip);
}
static const int32_t cp[]={0,2,4,6}, ri[]={0,1,1,2,0,2};
static const double v[]={4,2,5,3,1,6};
static const double dense[]={4,0,1,2,5,0,0,3,6};
static double complex zv[6];
static int factor(void *h,int z,int reuse) {
    if(z) return reuse ? fs_native_refactor_z(h,3,6,cp,ri,zv) : fs_native_factor_z(h,3,6,cp,ri,zv);
    return reuse ? fs_native_refactor_d(h,3,6,cp,ri,v) : fs_native_factor_d(h,3,6,cp,ri,v);
}
static void stats(void *h,struct fs_native_memory_stats *s) {CHECK(fs_native_get_memory_stats(h,s)==0);}
static void solve(void *h,int z,int op,int count,int expected,int alias) {
    double b[96],x[96];double complex zb[96],zx[96];
    for(int r=0;r<count;++r) for(int i=0;i<3;++i) {
        b[3*r+i]=0;zb[3*r+i]=0;
        for(int j=0;j<3;++j) {
            double a=op ? dense[j*3+i] : dense[i*3+j];
            double complex az=a*(1+0.2*I);if(op==2)az=conj(az);
            b[3*r+i]+=a*(j+r+1);
            zb[3*r+i]+=az*((j+r+1)+(j-1)*I);
        }
        x[3*r+i]=-91;zx[3*r+i]=-91+17*I;
    }
    double saved[96];double complex zsaved[96];
    memcpy(saved,b,sizeof saved);memcpy(zsaved,zb,sizeof zsaved);
    int code=z ? fs_native_solve_z(h,op,3,count,zb,alias?zb:zx) : fs_native_solve_d(h,op,3,count,b,alias?b:x);
    CHECK(code==expected);
    if(code==0) {
        for(int r=0;r<count;++r) for(int i=0;i<3;++i) {
            if(z) CHECK(cabs((alias?zb:zx)[3*r+i]-((i+r+1)+(i-1)*I))<1e-11);
            else CHECK(fabs((alias?b:x)[3*r+i]-(i+r+1))<1e-11);
        }
    } else {
        for(int r=0;r<count;++r) for(int i=0;i<3;++i) {
            if(z) CHECK((alias?zb:zx)[3*r+i]==(alias?zsaved[3*r+i]:-91+17*I));
            else CHECK((alias?b:x)[3*r+i]==(alias?saved[3*r+i]:-91));
        }
    }
}
static void one_kind(int z) {
    void *h=NULL,*other=NULL;
    struct fs_native_memory_stats s,base;
    int64_t calls,live,peak,count,ignored;
    const int64_t value_bytes=z?16:8, conversion=16+(4+value_bytes)*6;
    char message[160];
    CHECK(fs_native_create(z,&h)==0);
    stats(h,&s);CHECK(s.host_current_bytes==0&&s.host_peak_bytes==0&&s.memory_policy==0);
    CHECK(fs_native_set_memory(NULL,1,0,0)==2);
    CHECK(fs_native_get_memory_stats(h,NULL)==2);
    CHECK(fs_native_get_memory_stats(NULL,&s)==2);
    CHECK(factor(h,z,0)==0); /* Prepare the fixed backend object outside the payload oracle. */
    CHECK(fs_native_set_memory(h,1000000,7,FS_NATIVE_MEMORY_LOW)==0);
    fsparse_spd_probe_begin(0);
    CHECK(factor(h,z,0)==0);
    stats(h,&base);
    fsparse_spd_probe_end(&calls,&live,&peak);
    CHECK(live==base.host_current_bytes&&peak==base.host_peak_bytes);
    CHECK(base.strict_memory_supported==1&&base.peak_is_bound==1);
    CHECK(base.device_limit_bytes==7&&base.device_current_bytes==0&&base.device_peak_bytes==0);
    CHECK(base.host_peak_bytes>=base.host_current_bytes+conversion);
    count=calls;
    printf("C ABI %s payload current=%lld peak=%lld allocations=%lld\n",z?"complex":"real",(long long)live,(long long)peak,(long long)calls);
    /* Reject policy requests transactionally; querying stats preserves last_error. */
    CHECK(fs_native_set_memory(h,-1,0,0)==2);
    stats(h,&s);CHECK(s.host_limit_bytes==1000000);
    CHECK(fs_native_last_error(h,message,sizeof message)==2);
    solve(h,z,0,2,0,0);
    CHECK(fs_native_set_memory(h,1000000,7,9)==2);
    solve(h,z,1,2,0,1);
    CHECK(fs_native_set_memory(h,1000000,-1,0)==2);
    solve(h,z,2,2,0,0);
    CHECK(fs_native_set_memory(h,1000000,7,FS_NATIVE_MEMORY_LOW)==0);
    solve(h,z,0,2,0,0); /* Same policy retains factors. */
    /* All conversion/backend allocation failures cleanly invalidate and recover. */
    for(int64_t nth=1;nth<=count;++nth) {
        CHECK(fs_native_set_memory(h,1000000+2*nth,0,FS_NATIVE_MEMORY_LOW)==0);
        fsparse_spd_probe_begin(nth);
        CHECK(factor(h,z,0)==6);
        CHECK(fs_native_last_error(h,message,sizeof message)==6&&message[0]);
        CHECK(fs_native_set_memory(h,1000001+2*nth,0,FS_NATIVE_MEMORY_LOW)==0);
        fsparse_spd_probe_end(&ignored,&live,&peak);
        CHECK(live==0);
    }
    CHECK(factor(h,z,0)==0);
    solve(h,z,0,2,0,0);
    /* Tiny and exhausted conversion quotas never allocate or become unlimited. */
    for(int trial=0;trial<2;++trial) {
        CHECK(fs_native_set_memory(h,trial?conversion:1,0,FS_NATIVE_MEMORY_LOW)==0);
        fsparse_spd_probe_begin(0);
        CHECK(factor(h,z,0)==8);
        fsparse_spd_probe_end(&calls,&live,&peak);
        CHECK(calls==0&&live==0&&peak==0);
        stats(h,&s);CHECK(s.budget_required_bytes>s.budget_available_bytes);
        CHECK(fs_native_last_error(h,message,sizeof message)==8);
        solve(h,z,0,2,4,0);
    }
    /* The complete fresh LOW factor peak is a sufficient payload cap. */
    CHECK(fs_native_set_memory(h,base.host_peak_bytes-1,0,FS_NATIVE_MEMORY_LOW)==0);
    CHECK(factor(h,z,0)==8);
    stats(h,&s);CHECK(s.host_peak_bytes<=s.host_limit_bytes);
    CHECK(fs_native_set_memory(h,base.host_peak_bytes,0,FS_NATIVE_MEMORY_LOW)==0);
    CHECK(factor(h,z,0)==0);
    stats(h,&s);CHECK(s.host_current_bytes==base.host_current_bytes);
    /* Larger replacement conversion must coexist with the still-retained factors. */
    int32_t full_cp[]={0,3,6,9},full_ri[]={0,1,2,0,1,2,0,1,2};
    double full_v[]={4,2,0,0,5,3,1,0,6};double complex full_z[9];
    for(int i=0;i<9;++i)full_z[i]=full_v[i]*(1+0.2*I);
    fsparse_spd_probe_begin(0);
    int changed=z ? fs_native_refactor_z(h,3,9,full_cp,full_ri,full_z) : fs_native_refactor_d(h,3,9,full_cp,full_ri,full_v);
    fsparse_spd_probe_end(&calls,&live,&peak);
    CHECK(changed==8&&calls==0);
    stats(h,&s);CHECK(s.budget_required_bytes==16+(4+value_bytes)*9);
    CHECK(s.budget_available_bytes==base.host_peak_bytes-base.host_current_bytes);
    CHECK(fs_native_set_memory(h,1000000,0,FS_NATIVE_MEMORY_LOW)==0);
    changed=z ? fs_native_refactor_z(h,3,9,full_cp,full_ri,full_z) : fs_native_refactor_d(h,3,9,full_cp,full_ri,full_v);
    CHECK(changed==0);solve(h,z,0,2,0,0);
    /* Solve staging is transactional and accounted at its exact simultaneous peak. */
    int64_t solve_cap=base.host_current_bytes+value_bytes*3*32;
    CHECK(fs_native_set_memory(h,solve_cap,0,FS_NATIVE_MEMORY_LOW)==0);
    CHECK(factor(h,z,0)==0);
    for(int op=0;op<3;++op) {
        fsparse_spd_probe_begin(0);
        solve(h,z,op,32,0,op==2);
        fsparse_spd_probe_end(&calls,&live,&peak);
        CHECK(calls==1&&live==0&&peak==value_bytes*3*32);
        stats(h,&s);CHECK(s.host_peak_bytes==solve_cap&&s.peak_is_bound==0);
        CHECK(s.budget_required_bytes==0&&s.budget_available_bytes==0);
    }
    fsparse_spd_probe_begin(1);
    solve(h,z,0,32,6,1);
    fsparse_spd_probe_end(&calls,&live,&peak);CHECK(live==0);
    solve(h,z,0,1,0,0);
    CHECK(fs_native_set_memory(h,solve_cap-1,0,FS_NATIVE_MEMORY_LOW)==0);
    CHECK(factor(h,z,0)==0);
    solve(h,z,2,32,8,1);
    stats(h,&s);CHECK(s.budget_required_bytes==value_bytes*3*32);
    CHECK(s.budget_available_bytes==value_bytes*3*32-1);
    /* One failed solve does not destroy factors; smaller staging fits. */
    solve(h,z,1,1,0,0);
    CHECK(fs_native_create(z,&other)==0);
    CHECK(factor(other,z,0)==0);
    CHECK(fs_native_set_memory(h,1,0,FS_NATIVE_MEMORY_AUTO)==0);
    CHECK(factor(h,z,0)==8);
    solve(other,z,0,2,0,0);
    CHECK(fs_native_set_memory(h,1000000,0,FS_NATIVE_MEMORY_LOW)==0);
    fsparse_spd_probe_begin(0);
    CHECK(factor(h,z,0)==0);
    stats(h,&s);int64_t sequence_bound=s.host_peak_bytes;
    CHECK(factor(h,z,1)==0);
    stats(h,&s);
    fsparse_spd_probe_end(&calls,&live,&peak);
    if(s.host_peak_bytes>sequence_bound)sequence_bound=s.host_peak_bytes;
    CHECK(live==s.host_current_bytes&&peak<=sequence_bound);
    CHECK(s.host_peak_bytes<=s.host_limit_bytes);
    CHECK(fs_native_set_refine(h,1)==0);
    fsparse_spd_probe_begin(0);
    CHECK(factor(h,z,0)==0);
    stats(h,&s);fsparse_spd_probe_end(&calls,&live,&peak);
    CHECK(live==s.host_current_bytes&&peak==s.host_peak_bytes);
    printf("C ABI %s refined payload current=%lld peak=%lld allocations=%lld\n",z?"complex":"real",(long long)live,(long long)peak,(long long)calls);
    count=calls;
    for(int64_t nth=1;nth<=count;++nth) {
        CHECK(fs_native_set_memory(h,2000000+2*nth,0,FS_NATIVE_MEMORY_LOW)==0);
        fsparse_spd_probe_begin(nth);
        CHECK(factor(h,z,0)==6);
        CHECK(fs_native_set_memory(h,2000001+2*nth,0,FS_NATIVE_MEMORY_LOW)==0);
        fsparse_spd_probe_end(&ignored,&live,&peak);CHECK(live==0);
    }
    CHECK(factor(h,z,0)==0);
    for(int op=0;op<3;++op)solve(h,z,op,2,0,op==1);
    CHECK(fs_native_set_memory(h,0,0,FS_NATIVE_MEMORY_AUTO)==0);
    CHECK(factor(h,z,1)==0);solve(h,z,0,2,0,0);
    if(fs_native_ordering_available(FS_NATIVE_ORDER_AMD)==0) {
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_AMD)==0);
        CHECK(factor(h,z,0)==0);stats(h,&s);
        CHECK(s.host_peak_bytes==-1&&s.strict_memory_supported==0&&s.peak_is_bound==1);
        CHECK(fs_native_set_memory(h,1000000,0,FS_NATIVE_MEMORY_LOW)==0);
        int result=factor(h,z,0);
        CHECK(result==(fsparse_native_amd_bounded_available()?0:3));
        if(result==0) {
            stats(h,&s);CHECK(s.strict_memory_supported==1&&s.host_peak_bytes<=s.host_limit_bytes);
            for(int op=0;op<3;++op)solve(h,z,op,2,0,0);
        }
    }
    /* ND is actual optional graph-library execution; its opaque storage is uncapped only. */
    CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_DEFAULT)==0);
    CHECK(fs_native_set_memory(h,0,0,FS_NATIVE_MEMORY_AUTO)==0);
    CHECK(factor(h,z,0)==0);
    int nd=fs_native_ordering_available(FS_NATIVE_ORDER_ND), before_nd=nd_calls;
    CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_ND)==nd);
    if(nd==0) {
        solve(h,z,0,2,4,0); /* Successful ordering change invalidates. */
        CHECK(factor(h,z,0)==0&&nd_calls>before_nd);
        stats(h,&s);CHECK(s.host_peak_bytes==-1&&s.strict_memory_supported==0&&s.peak_is_bound==1);
        for(int op=0;op<3;++op)solve(h,z,op,2,0,op==1);
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_ND)==0);
        solve(h,z,0,2,0,0); /* Same order retains actual ND factors. */
        CHECK(factor(h,z,1)==0);solve(h,z,2,2,0,0);
        before_nd=nd_calls;
        CHECK(fs_native_set_memory(h,1000000,0,FS_NATIVE_MEMORY_AUTO)==0);
        CHECK(factor(h,z,0)==3&&nd_calls==before_nd);
        stats(h,&s);CHECK(s.strict_memory_supported==0);
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_DEFAULT)==0);
        CHECK(factor(h,z,0)==0);solve(h,z,0,2,0,0);
        CHECK(fs_native_set_memory(h,0,0,FS_NATIVE_MEMORY_LOW)==0);
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_ND)==0);
        CHECK(factor(h,z,0)==3&&nd_calls==before_nd);
        stats(h,&s);CHECK(s.strict_memory_supported==0);
        CHECK(fs_native_set_ordering(h,FS_NATIVE_ORDER_DEFAULT)==0);
        CHECK(factor(h,z,0)==0);solve(h,z,1,2,0,0);
    } else {
        CHECK(nd==3&&nd_calls==before_nd);
        solve(h,z,0,2,0,0); /* Unavailable setter retains default factors. */
    }
    CHECK(fs_native_free(&h)==0&&fs_native_free(&h)==0&&h==NULL);
    CHECK(fs_native_free(&other)==0&&other==NULL);
}
int main(void) {
#ifdef FS_TEST_EXPECT_METIS
    CHECK(fs_native_ordering_available(FS_NATIVE_ORDER_ND)==0);
#else
    CHECK(fs_native_ordering_available(FS_NATIVE_ORDER_ND)==3);
#endif
    for(int i=0;i<6;++i)zv[i]=v[i]*(1+0.2*I);
    one_kind(0);one_kind(1);
    puts("C ABI complete payload caps, staging, allocation failure and independent handles pass");
}
