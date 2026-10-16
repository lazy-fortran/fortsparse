#include "fsparse_spd_cuda.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace {
void require(bool ok) { if (!ok) { std::fputs("CUDA DAG oracle failed\n",stderr); std::exit(1); } }
void check(int code) { if (code) { std::fprintf(stderr,"CUDA DAG status %d\n",code); std::exit(1); } }
void close(double a, double b) { require(std::isfinite(a) && std::abs(a-b)<3e-13*std::max(1.,std::abs(b))); }
constexpr int n = 6;
constexpr size_t capacity = 64;
const fsparse_spd_cuda_level_panel panels[] = {
    {3,1,5,7,0}, {12,1,5,6,0}, {20,1,4,6,1}, {28,3,3,5,2}};
const int global_rows[4][5] = {{0,2,3,4,5},{1,2,3,4,5},{2,3,4,5},{3,4,5}};
const int first_columns[] = {0,1,2,3};
const int maps[] = {1,2,3,4, 1, 1,2,3, 1,2,3};
fsparse_spd_cuda_level_task make_task(size_t source, size_t dest, int first,
                                     int m, int columns, size_t row_map, size_t col_map)
{
    const auto &s = panels[source], &d = panels[dest];
    return {{s.offset+first,s.offset+first,d.offset,row_map,col_map,
             s.lda,s.lda,m,columns,s.width,d.lda,d.width},source,dest};
}
std::vector<fsparse_spd_cuda_level_task> tasks()
{
    return {make_task(0,2,1,4,1,0,4), make_task(0,3,2,3,3,5,8),
            make_task(1,2,1,4,1,0,4), make_task(1,3,2,3,3,5,8),
            make_task(2,3,1,3,3,5,8)};
}

void manufactured(int repeat, std::vector<double> &input,
                  std::vector<double> &lower, std::vector<double> &matrix)
{
    lower.assign(n*n,0); matrix.assign(n*n,0); input.assign(capacity,9917);
    for (int j=0;j<n;++j)
        for (int i=j;i<n;++i)
            lower[i+n*j] = i==j ? 2+.125*j : .03125*(1+(i+2*j+repeat)%7);
    lower[1] = 0; // Independent first two sources, both feed all later panels.
    for (int j=0;j<n;++j)
        for (int i=0;i<n;++i)
            for (int k=0;k<n;++k)
                matrix[i+n*j] += lower[i+n*k]*lower[j+n*k];
    for (int p=0;p<4;++p)
        for (int j=0;j<panels[p].width;++j)
            for (int i=j;i<panels[p].rows;++i)
                input[panels[p].offset+i+j*panels[p].lda] =
                    matrix[global_rows[p][i]+n*(first_columns[p]+j)];
}

void run_levels(void *ctx, int &info)
{
    check(fsparse_spd_cuda_begin_factor(ctx));
    for (size_t level=0;level<3;++level) {
        check(fsparse_spd_cuda_factor_level(ctx,level));
        check(fsparse_spd_cuda_update_level(ctx,level));
    }
    check(fsparse_spd_cuda_factor_info(ctx,&info));
}

void independent_oracle(const std::vector<double> &actual,
                         const std::vector<double> &input,
                         const std::vector<double> &expected,
                         const std::vector<double> &matrix)
{
    std::vector<double> reconstructed(n*n,0);
    std::vector<bool> changed(capacity,false);
    for (int p=0;p<4;++p)
        for (int j=0;j<panels[p].width;++j)
            for (int i=j;i<panels[p].rows;++i) {
                const size_t entry=panels[p].offset+i+j*panels[p].lda;
                const int row=global_rows[p][i], col=first_columns[p]+j;
                close(actual[entry],expected[row+n*col]);
                reconstructed[row+n*col]=actual[entry]; changed[entry]=true;
            }
    for (size_t i=0;i<capacity;++i) if (!changed[i]) close(actual[i],input[i]);
    // Dense independent solve and original-A backward error, with known x.
    double x[n], b[n]={}, solution[n], temporary[n];
    for (int i=0;i<n;++i) x[i]=i%2 ? -.5 : 1+i*.125;
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) b[i]+=matrix[i+n*j]*x[j];
        double value=b[i];
        for (int j=0;j<i;++j) value-=reconstructed[i+n*j]*temporary[j];
        temporary[i]=value/reconstructed[i+n*i];
    }
    for (int i=n-1;i>=0;--i) {
        double value=temporary[i];
        for (int j=i+1;j<n;++j) value-=reconstructed[j+n*i]*solution[j];
        solution[i]=value/reconstructed[i+n*i]; close(solution[i],x[i]);
    }
    double residual=0, anorm=0, xnorm=0, bnorm=0;
    for (int i=0;i<n;++i) {
        double row_sum=0, value=-b[i];
        for (int j=0;j<n;++j) {
            row_sum+=std::abs(matrix[i+n*j]); value+=matrix[i+n*j]*solution[j];
        }
        residual=std::max(residual,std::abs(value)); anorm=std::max(anorm,row_sum);
        xnorm=std::max(xnorm,std::abs(solution[i])); bnorm=std::max(bnorm,std::abs(b[i]));
    }
    require(residual/(anorm*xnorm+bnorm)<1e-12);
}

void invalid_plans(void *ctx, const std::vector<fsparse_spd_cuda_level_task> &good)
{
    auto t=good;
    auto p=std::vector<fsparse_spd_cuda_level_panel>(panels,panels+4);
    t[0].destination=0;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    t=good; p[2].level=0;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    p.assign(panels,panels+4); p[1].offset=4;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    p.assign(panels,panels+4); t[0].source=4;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    t=good; t[0].task.left_offset=capacity;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    t=good; t[0].task.row_map_offset=100;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    t=good; t[0].task.destination_rows=4;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    p[0].lda=0;
    require(fsparse_spd_cuda_set_levels(ctx,p.data(),4,t.data(),t.size(),3)==-1);
    require(fsparse_spd_cuda_set_levels(ctx,nullptr,1,nullptr,0,1)==-1);
    require(fsparse_spd_cuda_factor_level(ctx,3)==-1);
    require(fsparse_spd_cuda_update_level(ctx,3)==-1);
}
}
int main()
{
    void *ctx=nullptr;
    check(fsparse_spd_cuda_create_deferred(&ctx,0,capacity));
    check(fsparse_spd_cuda_set_levels(ctx,nullptr,0,nullptr,0,0));
    check(fsparse_spd_cuda_set_maps(ctx,maps,sizeof(maps)/sizeof(*maps)));
    const size_t before=fsparse_spd_cuda_allocated_bytes(ctx);
    auto plan=tasks();
    check(fsparse_spd_cuda_set_levels(ctx,panels,4,plan.data(),plan.size(),3));
    const size_t retained=fsparse_spd_cuda_allocated_bytes(ctx);
    const size_t host_retained=fsparse_spd_cuda_host_plan_bytes(ctx);
    require(retained>before && host_retained>0);
    invalid_plans(ctx,plan); // Rejected setup must preserve the last valid plan.
    std::vector<double> input,expected,matrix,actual(capacity);
    for (int repeat=0;repeat<3;++repeat) {
        manufactured(repeat,input,expected,matrix);
        if (repeat==1) {
            auto bad=input; bad[panels[0].offset]=-1;
            check(fsparse_spd_cuda_upload(ctx,0,bad.data(),capacity));
            int info=0; run_levels(ctx,info); require(info!=0);
            // Other same-level panels may have executed. All failed output is
            // discarded; the next whole refactor must clear failure state.
        }
        check(fsparse_spd_cuda_upload(ctx,0,input.data(),capacity));
        int info=-1; run_levels(ctx,info); require(info==0);
        check(fsparse_spd_cuda_download(ctx,0,actual.data(),capacity));
        check(fsparse_spd_cuda_synchronize(ctx));
        independent_oracle(actual,input,expected,matrix);
        require(fsparse_spd_cuda_allocated_bytes(ctx)==retained);
        require(fsparse_spd_cuda_host_plan_bytes(ctx)==host_retained);
    }
    // Updating maps invalidates retained schedules rather than using stale IDs.
    check(fsparse_spd_cuda_set_maps(ctx,maps,sizeof(maps)/sizeof(*maps)));
    require(fsparse_spd_cuda_allocated_bytes(ctx)==before);
    check(fsparse_spd_cuda_begin_factor(ctx));
    require(fsparse_spd_cuda_factor_level(ctx,0)==-1);
    int info=0; check(fsparse_spd_cuda_factor_info(ctx,&info));
    check(fsparse_spd_cuda_destroy(ctx));
    std::puts("CUDA DAG shared destinations, chains, refactor, failure and bounds pass");
}
