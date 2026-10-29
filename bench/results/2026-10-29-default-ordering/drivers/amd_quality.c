#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <amd.h>
#include "fsparse_native_amd_bounded.h"
#include "fsparse_native_amd_core.h"
static long fill(int n,const int*ptr,const int*row,const int*perm,const int*inv){
  int *parent=malloc(n*sizeof(int)),*mark=malloc(n*sizeof(int)); long cnt=0;
  for(int k=0;k<n;k++){parent[k]=-1;mark[k]=-1;}
  for(int k=0;k<n;k++){int old=perm[k]-1; mark[k]=k; cnt++;
    for(int p=ptr[old]-1;p<ptr[old+1]-1;p++){int i=inv[row[p]-1]-1; if(i>=k)continue;
      while(mark[i]!=k){mark[i]=k;cnt++; if(parent[i]==-1){parent[i]=k;break;} i=parent[i];}}}
  free(parent);free(mark);return cnt;}
static unsigned long long rs=88172645463325252ULL; static unsigned rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return (unsigned)(rs>>11);}
/* symmetric CSC (1-based, with diagonal) from edge list */
static void build(int n,int ne,int*ei,int*ej,int**ptr,int**row,int*nnz){
  char*seen=NULL; int *deg=calloc(n,sizeof(int));
  /* use sort-based dedupe: collect pairs */
  long m=2L*ne+n; int64_t*key=malloc(m*sizeof(int64_t)); long q=0;
  for(int e=0;e<ne;e++){ if(ei[e]==ej[e])continue; key[q++]=(int64_t)ej[e]*n+ei[e]; key[q++]=(int64_t)ei[e]*n+ej[e];}
  for(int i=0;i<n;i++)key[q++]=(int64_t)i*n+i;
  int cmp(const void*a,const void*b){int64_t x=*(int64_t*)a,y=*(int64_t*)b;return x<y?-1:x>y;}
  qsort(key,q,sizeof(int64_t),cmp); long u=0; for(long t=0;t<q;t++) if(!u||key[t]!=key[u-1])key[u++]=key[t];
  *ptr=calloc(n+1,sizeof(int)); *row=malloc(u*sizeof(int));
  for(long t=0;t<u;t++){int c=key[t]/n; (*ptr)[c+1]++; (*row)[t]=key[t]%n+1;}
  (*ptr)[0]=1; for(int c=0;c<n;c++)(*ptr)[c+1]+=(*ptr)[c]; *nnz=u; free(key);free(deg);(void)seen;
}
static int run(const char*name,int n,int nnz,int*ptr,int*row){
  int *perm=malloc(n*sizeof(int)),*inv=malloc(n*sizeof(int));
  for(int i=0;i<n;i++){perm[i]=i+1;inv[i]=i+1;}
  long nat=fill(n,ptr,row,perm,inv); long own[2];
  for(int low=0;low<2;low++){int64_t ab,ob; int c=fsparse_native_amd_bounded_size(n,nnz,ptr,row,low,&ab,&ob);
    int32_t*arena=malloc(ab); c=c?c:fsparse_native_amd_bounded(n,nnz,ptr,row,low,arena,ab,perm,inv); free(arena);
    if(c){printf("%s low=%d FAILED code %d\n",name,low,c);return 1;}
    own[low]=fill(n,ptr,row,perm,inv);}
  int *ap=malloc((n+1)*sizeof(int)),*ai=malloc(nnz*sizeof(int)),*p=malloc(n*sizeof(int));
  for(int i=0;i<=n;i++)ap[i]=ptr[i]-1; for(int i=0;i<nnz;i++)ai[i]=row[i]-1;
  amd_order(n,ap,ai,p,NULL,NULL); for(int k=0;k<n;k++){perm[k]=p[k]+1;inv[p[k]]=k+1;}
  long ss=fill(n,ptr,row,perm,inv);
  printf("%-22s n=%7d nnz=%8d natural=%11ld own=%10ld own_low=%10ld suitesparse=%10ld own/ss=%.3f\n",name,n,nnz,nat,own[0],own[1],ss,(double)own[0]/ss);
  free(ap);free(ai);free(p);free(perm);free(inv); return own[0]!=own[1];
}
int main(void){
  int bad=0;
  /* grid 5-point, random relabelled */
  for(int s=0;s<3;s++){int nx=(int[]){50,200,400}[s]; int n=nx*nx, ne=0; int*ei=malloc(4*n*sizeof(int)),*ej=malloc(4*n*sizeof(int)); int*lab=malloc(n*sizeof(int));
    for(int i=0;i<n;i++)lab[i]=i; for(int i=n-1;i>0;i--){int j=rnd()%(i+1);int t=lab[i];lab[i]=lab[j];lab[j]=t;}
    for(int y=0;y<nx;y++)for(int x=0;x<nx;x++){int c=y*nx+x; if(x+1<nx){ei[ne]=lab[c];ej[ne++]=lab[c+1];} if(y+1<nx){ei[ne]=lab[c];ej[ne++]=lab[c+nx];}}
    int*ptr,*row,nnz; build(n,ne,ei,ej,&ptr,&row,&nnz); char nm[64]; sprintf(nm,"grid%d-shuffled",nx); bad|=run(nm,n,nnz,ptr,row);}
  /* 3D grid */
  {int nx=40,n=nx*nx*nx,ne=0; int*ei=malloc(3*n*sizeof(int)),*ej=malloc(3*n*sizeof(int));
   for(int z=0;z<nx;z++)for(int y=0;y<nx;y++)for(int x=0;x<nx;x++){int c=(z*nx+y)*nx+x; if(x+1<nx){ei[ne]=c;ej[ne++]=c+1;} if(y+1<nx){ei[ne]=c;ej[ne++]=c+nx;} if(z+1<nx){ei[ne]=c;ej[ne++]=c+nx*nx;}}
   int*ptr,*row,nnz; build(n,ne,ei,ej,&ptr,&row,&nnz); bad|=run("grid3d-40",n,nnz,ptr,row);}
  /* random sparse graphs, several densities, with some dense rows and isolated nodes */
  for(int t=0;t<40;t++){int n=1+rnd()%3000; int ne=rnd()%(6*n+1); int extra=(t%4==0)?3:0; int tot=ne+extra*n;
    int*ei=malloc((tot+1)*sizeof(int)),*ej=malloc((tot+1)*sizeof(int)); int k=0;
    for(int e=0;e<ne;e++){ei[k]=rnd()%n;ej[k++]=rnd()%n;}
    for(int d=0;d<extra;d++){int h=rnd()%n; for(int i=0;i<n;i++){ei[k]=h;ej[k++]=i;}}
    int*ptr,*row,nnz; build(n,k,ei,ej,&ptr,&row,&nnz); char nm[64]; sprintf(nm,"random%d",t);
    int*perm=malloc(n*sizeof(int)),*inv=malloc(n*sizeof(int));
    for(int low=0;low<2;low++){int64_t ab,ob; int c=fsparse_native_amd_bounded_size(n,nnz,ptr,row,low,&ab,&ob); int32_t*arena=malloc(ab);
      c=c?c:fsparse_native_amd_bounded(n,nnz,ptr,row,low,arena,ab,perm,inv); free(arena); if(c){printf("%s low=%d n=%d FAILED %d\n",nm,low,n,c);bad=1;}}
    if(t%8==0) bad|=run(nm,n,nnz,ptr,row);
  }
  /* one dense column and a chain */
  printf(bad?"STRESS FAIL\n":"STRESS PASS\n"); return bad;
}
