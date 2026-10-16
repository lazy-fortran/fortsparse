#include "fsparse_spd_cuda.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>
namespace {
void require(bool ok) { if (!ok) { std::fputs("GPU triangular solve oracle failed\n",stderr); std::exit(1); } }
void check(int code) { if (code) { std::fprintf(stderr,"GPU solve code %d\n",code); std::exit(1); } }
void close(double a,double b) { require(std::isfinite(a) && std::abs(a-b)<4e-13*std::max(1.,std::abs(b))); }

void verify(void *ctx, int n, const std::vector<double> &lower,
            const std::vector<double> &input, size_t rhs_offset)
{
    const int ld=n+3;
    std::vector<double> matrix(n*n,0), b(ld*8,9937), actual(ld*8), expected(n*8);
    for (int j=0;j<n;++j)
        for (int i=0;i<n;++i)
            for (int k=0;k<n;++k) matrix[i+n*j]+=lower[i+n*k]*lower[j+n*k];
    for (int r=0;r<8;++r) {
        for (int i=0;i<n;++i) expected[i+n*r]=(r==4)?0.0:((i+2*r)%13-6)*.125;
        for (int i=0;i<n;++i) {
            b[i+ld*r]=0;
            for (int j=0;j<n;++j) b[i+ld*r]+=matrix[i+n*j]*expected[j+n*r];
        }
    }
    const size_t device=fsparse_spd_cuda_allocated_bytes(ctx),host=fsparse_spd_cuda_host_plan_bytes(ctx);
    for (int nrhs : {0,1,8}) {
        for (int repeat=0;repeat<2;++repeat) {
            check(fsparse_spd_cuda_upload(ctx,rhs_offset,b.data(),b.size()));
            check(fsparse_spd_cuda_solve(ctx,rhs_offset,ld,nrhs));
            check(fsparse_spd_cuda_download(ctx,rhs_offset,actual.data(),actual.size()));
            check(fsparse_spd_cuda_synchronize(ctx));
            for (int r=0;r<8;++r)
                for (int i=0;i<ld;++i) close(actual[i+ld*r],r<nrhs && i<n?expected[i+n*r]:b[i+ld*r]);
            for (int r=0;r<nrhs;++r) {
                double residual=0,anorm=0,xnorm=0,bnorm=0;
                for (int i=0;i<n;++i) {
                    double value=-b[i+ld*r],row_sum=0;
                    for (int j=0;j<n;++j) { value+=matrix[i+n*j]*actual[j+ld*r]; row_sum+=std::abs(matrix[i+n*j]); }
                    residual=std::max(residual,std::abs(value)); anorm=std::max(anorm,row_sum);
                    xnorm=std::max(xnorm,std::abs(actual[i+ld*r])); bnorm=std::max(bnorm,std::abs(b[i+ld*r]));
                }
                require(residual<=1e-12*std::max(1.,anorm*xnorm+bnorm));
            }
        }
    }
    std::vector<double> unchanged(input.size());
    check(fsparse_spd_cuda_download(ctx,0,unchanged.data(),unchanged.size()));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (size_t i=0;i<input.size();++i)
        if (std::isnan(input[i])) require(std::isnan(unchanged[i])); else close(unchanged[i],input[i]);
    require(device==fsparse_spd_cuda_allocated_bytes(ctx) && host==fsparse_spd_cuda_host_plan_bytes(ctx));
    for (int code : {fsparse_spd_cuda_solve(nullptr,rhs_offset,ld,1),
        fsparse_spd_cuda_solve(ctx,rhs_offset,n-1,1), fsparse_spd_cuda_solve(ctx,rhs_offset,ld,-1),
        fsparse_spd_cuda_solve(ctx,rhs_offset,ld,9),fsparse_spd_cuda_solve(ctx,0,ld,1),
        fsparse_spd_cuda_solve(ctx,std::numeric_limits<size_t>::max(),ld,1)}) require(code==-1);
}

void branch_case()
{
    constexpr int n=6;
    const fsparse_spd_cuda_level_panel panels[]={{3,1,5,7,0},{12,1,5,6,0},{20,1,4,6,1},{28,3,3,5,2}};
    const int tails[]={3,4,5,6,3,4,5,6,4,5,6};
    const int rows[4][5]={{0,2,3,4,5},{1,2,3,4,5},{2,3,4,5},{3,4,5}};
    const int columns[]={0,1,2,3};
    constexpr size_t rhs_offset=64;
    void *ctx=nullptr;
    check(fsparse_spd_cuda_create_deferred(&ctx,0,256));
    check(fsparse_spd_cuda_set_solve_plan(ctx,n,panels,4,tails,11,3));
    require(fsparse_spd_cuda_solve(ctx,rhs_offset,n+3,1)==-1); // no completed factor
    for (int pass=0;pass<3;++pass) {
        std::vector<double> lower(n*n,0),input(rhs_offset,9937);
        for (int j=0;j<n;++j)
            for (int i=j;i<n;++i) lower[i+n*j]=i==j?2+.125*(j+pass):.03125*(1+(i+2*j+pass)%7);
        lower[1]=0; // same-level sources with overlapping future destinations
        for (int p=0;p<4;++p)
            for (int j=0;j<panels[p].width;++j)
                for (int i=0;i<panels[p].rows;++i)
                    input[panels[p].offset+i+j*panels[p].lda]=i<j?
                        std::numeric_limits<double>::quiet_NaN():lower[rows[p][i]+n*(columns[p]+j)];
        // An independent producer supplies known L. factor_info finalizes its
        // caller-managed factor stage; the production supernode test factors A.
        check(fsparse_spd_cuda_begin_factor(ctx));
        check(fsparse_spd_cuda_upload(ctx,0,input.data(),input.size()));
        int info=-1; check(fsparse_spd_cuda_factor_info(ctx,&info)); require(info==0);
        verify(ctx,n,lower,input,rhs_offset);
    }
    // Invalid metadata must not replace the valid plan.
    int invalid[11]; std::copy(tails,tails+11,invalid); invalid[0]=0;
    require(fsparse_spd_cuda_set_solve_plan(ctx,n,panels,4,invalid,11,3)==-1);
    std::copy(tails,tails+11,invalid); invalid[1]=3;
    require(fsparse_spd_cuda_set_solve_plan(ctx,n,panels,4,invalid,11,3)==-1);
    auto wrong=std::vector<fsparse_spd_cuda_level_panel>(panels,panels+4); wrong[2].level=0;
    require(fsparse_spd_cuda_set_solve_plan(ctx,n,wrong.data(),4,tails,11,3)==-1);
    require(fsparse_spd_cuda_set_solve_plan(nullptr,n,panels,4,tails,11,3)==-1);
    require(fsparse_spd_cuda_set_solve_plan(ctx,n,panels,4,tails,10,3)==-1);
    wrong.assign(panels,panels+4); wrong[1].offset=3;
    require(fsparse_spd_cuda_set_solve_plan(ctx,n,wrong.data(),4,tails,11,3)==-1);
    check(fsparse_spd_cuda_solve(ctx,rhs_offset,n+3,0));
    const double bad=-1;
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_upload(ctx,3,&bad,1));
    check(fsparse_spd_cuda_panel_factor(ctx,3,1,7));
    int info=0; check(fsparse_spd_cuda_factor_info(ctx,&info)); require(info==1);
    require(fsparse_spd_cuda_solve(ctx,rhs_offset,n+3,1)==-1);
    check(fsparse_spd_cuda_set_maps(ctx,nullptr,0));
    require(fsparse_spd_cuda_solve(ctx,rhs_offset,n+3,1)==-1);
    check(fsparse_spd_cuda_destroy(ctx));
}

void dense_case(int width)
{
    const int ld=width+2;
    const size_t rhs_offset=3+size_t(ld)*width;
    const fsparse_spd_cuda_level_panel panel={3,width,width,ld,0};
    std::vector<double> lower(width*width,0),input(rhs_offset,9937);
    for (int j=0;j<width;++j)
        for (int i=0;i<width;++i) {
            lower[i+width*j]=i==j?2+.001*j:(i>j?.001*((i+j)%7-3):0);
            input[3+i+size_t(j)*ld]=i>=j?lower[i+width*j]:std::numeric_limits<double>::quiet_NaN();
        }
    void *ctx=nullptr;
    check(fsparse_spd_cuda_create_deferred(&ctx,0,rhs_offset+size_t(width+3)*8));
    check(fsparse_spd_cuda_set_solve_plan(ctx,width,&panel,1,nullptr,0,1));
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_upload(ctx,0,input.data(),input.size()));
    int info=-1; check(fsparse_spd_cuda_factor_info(ctx,&info)); require(info==0);
    verify(ctx,width,lower,input,rhs_offset);
    check(fsparse_spd_cuda_destroy(ctx));
}
}
int main()
{
    branch_case();
    for (int width : {1,8,32,128}) dense_case(width);
    std::puts("GPU level triangular solves, RHS tiles, padding, reuse and failures pass");
}
