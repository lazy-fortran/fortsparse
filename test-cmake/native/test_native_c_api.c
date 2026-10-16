#include "fortsparse_native.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

/* Independent original dense A actions generate RHS for known x, never using
 * the CSC converter or solver. Two RHS distinguish column-major strides. */
static void real_rhs(const double *a, int op, double *b) {
    for (int r=0;r<2;++r) for (int i=0;i<3;++i) {
        b[r*3+i]=0;
        for(int j=0;j<3;++j) b[r*3+i]+=(op ? a[j*3+i] : a[i*3+j])*(j+1+r);
    }
}
static void complex_rhs(const double complex *a, int op, double complex *b) {
    for (int r=0;r<2;++r) for (int i=0;i<3;++i) {
        b[r*3+i]=0;
        for(int j=0;j<3;++j) {
            double complex v=op ? a[j*3+i] : a[i*3+j];
            if(op==2) v=conj(v);
            b[r*3+i]+=v*((j+1+r)+(2-j)*I);
        }
    }
}
int main(void) {
    void *d=NULL,*z=NULL;
    int32_t cp[]={0,2,4,6},ri[]={0,1,1,2,0,2};
    double v[]={4,2,5,3,1,6};
    double a[]={4,0,1,2,5,0,0,3,6},b[6],x[6];
    double complex vz[]={4+I,2-I,5+2*I,3+I,1-I,6-I};
    double complex az[]={4+I,0,1-I,2-I,5+2*I,0,0,3+I,6-I},bz[6],xz[6];
    char msg[8];
    CHECK(fs_native_solve_d(NULL,0,3,1,b,x)==2);
    CHECK(fs_native_create(0,NULL)==2);
    CHECK(fs_native_free(NULL)==0);
    CHECK(fs_native_create(7,&d)==2 && d==NULL);
    CHECK(fs_native_create(0,&d)==0 && d!=NULL);
    CHECK(fs_native_create(0,&d)==2); /* do not leak/overwrite live handle */
    CHECK(fs_native_create(1,&z)==0);
    CHECK(fs_native_solve_d(d,0,3,1,b,x)==4);
    CHECK(fs_native_set_refine(d,1)==0);
    CHECK(fs_native_set_refine(d,3)==2);
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==0);
    CHECK(fs_native_factor_z(z,3,6,cp,ri,vz)==0);
    v[0]=999; vz[0]=999; /* borrowed input storage can change after factor */
    for (int op=0;op<3;++op) {
        real_rhs(a,op,b); complex_rhs(az,op,bz);
        CHECK(fs_native_solve_d(d,op,3,2,b,x)==0);
        CHECK(fs_native_solve_z(z,op,3,2,bz,xz)==0);
        for(int r=0;r<2;++r) for(int i=0;i<3;++i) {
            CHECK(fabs(x[r*3+i]-(i+1+r))<1e-12);
            CHECK(cabs(xz[r*3+i]-((i+1+r)+(2-i)*I))<1e-12);
        }
        CHECK(fs_native_solve_d(d,op,3,2,b,b)==0); /* supported alias */
        CHECK(fs_native_solve_z(z,op,3,2,bz,bz)==0);
    }
    x[0]=123; b[0]=NAN;
    CHECK(fs_native_solve_d(d,0,3,2,b,x)==2 && x[0]==123);
    CHECK(fs_native_solve_d(d,8,3,2,b,x)==2 && x[0]==123);
    CHECK(fs_native_solve_d(d,0,2,2,b,x)==2);
    CHECK(fs_native_solve_d(d,0,3,INT_MAX,b,x)==7);
    CHECK(fs_native_last_error(d,msg,sizeof msg)==7 && msg[7]=='\0');
    CHECK(fs_native_solve_d(d,0,3,1,NULL,x)==2);
    v[0]=8; a[0]=8; vz[0]=8+I; az[0]=8+I;
    CHECK(fs_native_refactor_d(d,3,6,cp,ri,v)==0);
    CHECK(fs_native_refactor_z(z,3,6,cp,ri,vz)==0);
    real_rhs(a,0,b);complex_rhs(az,2,bz);
    CHECK(fs_native_solve_d(d,0,3,2,b,x)==0);
    CHECK(fs_native_solve_z(z,2,3,2,bz,xz)==0);
    CHECK(fabs(x[0]-1)<1e-12 && cabs(xz[0]-(1+2*I))<1e-12);
    ri[0]=3;
    CHECK(fs_native_refactor_d(d,3,6,cp,ri,v)==2);
    CHECK(fs_native_solve_d(d,0,3,2,b,x)==4);
    ri[0]=0;
    CHECK(fs_native_factor_d(d,INT_MAX,6,cp,ri,v)==7);
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==0);
    CHECK(fs_native_factor_z(d,3,6,cp,ri,vz)==2);
    CHECK(fs_native_solve_d(d,0,3,2,b,x)==4);
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==0);
    v[0]=INFINITY;
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==2);
    v[0]=8;
    ri[1]=0; /* duplicated row */
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==2);ri[1]=1;
    cp[1]=7;
    CHECK(fs_native_factor_d(d,3,6,cp,ri,v)==2);cp[1]=2;
    /* New pattern, singular refactor, recovery to unrelated diagonal shape. */
    int32_t dc[]={0,1,2,3},dr[]={0,1,2};double dv[]={1,0,3};
    CHECK(fs_native_refactor_d(d,3,3,dc,dr,dv)==1);
    CHECK(fs_native_solve_d(d,0,3,2,b,x)==4);
    dv[1]=2;
    CHECK(fs_native_refactor_d(d,3,3,dc,dr,dv)==0);
    b[0]=1;b[1]=4;b[2]=9;
    CHECK(fs_native_solve_d(d,0,3,1,b,x)==0 && fabs(x[2]-3)<1e-14);
    CHECK(fs_native_set_refine(d,0)==0);
    CHECK(fs_native_solve_d(d,0,3,1,b,x)==4);
    CHECK(fs_native_free(&d)==0 && d==NULL);
    CHECK(fs_native_free(&d)==0);
    CHECK(fs_native_free(&z)==0 && z==NULL);
    for(int cycle=0;cycle<50;++cycle) {
        CHECK(fs_native_create(0,&d)==0);
        CHECK(fs_native_factor_d(d,3,3,dc,dr,dv)==0);
        CHECK(fs_native_free(&d)==0);
    }
    return 0;
}
