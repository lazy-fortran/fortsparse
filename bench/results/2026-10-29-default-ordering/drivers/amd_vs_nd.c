#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <metis.h>
#include <amd.h>
#include "fsparse_native_amd_bounded.h"
/* Cholesky column counts of the symmetric pattern under perm: nnz(L) and flops sum c^2 */
static void stats(int n,const int*ptr,const int*row,const int*perm,const int*inv,long*nz,double*fl){
  int *parent=malloc(n*sizeof(int)),*mark=malloc(n*sizeof(int)); long *cc=calloc(n,sizeof(long)); long cnt=0;
  for(int k=0;k<n;k++){parent[k]=-1;mark[k]=-1;}
  for(int k=0;k<n;k++){int old=perm[k]-1; mark[k]=k; cnt++; cc[k]++;
    for(int p=ptr[old]-1;p<ptr[old+1]-1;p++){int i=inv[row[p]-1]-1; if(i>=k)continue;
      while(mark[i]!=k){mark[i]=k;cnt++;cc[i]++; if(parent[i]==-1){parent[i]=k;break;} i=parent[i];}}}
  double f=0; for(int k=0;k<n;k++) f+=(double)cc[k]*cc[k];
  *nz=cnt;*fl=f; free(parent);free(mark);free(cc);}
int main(int argc,char**argv){
  int m=atoi(argv[1]); int nv=m*m,n=2*nv;
  /* vertex adjacency of the P1 triangulation, 2 dofs/vertex -> dof pattern */
  int dx[7]={0,-1,1,0,0,1,-1},dy[7]={0,0,0,-1,1,-1,1};
  int*ptr=malloc((n+1)*sizeof(int)),*row=malloc((size_t)n*14*sizeof(int)); int q=0;
  for(int y=0;y<m;y++)for(int x=0;x<m;x++)for(int a=0;a<2;a++){int c=2*(y*m+x)+a; ptr[c]=q+1; int buf[14],nb=0;
    for(int t=0;t<7;t++){int xx=x+dx[t],yy=y+dy[t]; if(xx<0||yy<0||xx>=m||yy>=m)continue; for(int b=0;b<2;b++)buf[nb++]=2*(yy*m+xx)+b+1;}
    for(int i=1;i<nb;i++){int v=buf[i],j=i-1;while(j>=0&&buf[j]>v){buf[j+1]=buf[j];j--;}buf[j+1]=v;}
    for(int i=0;i<nb;i++)row[q++]=buf[i];}
  ptr[n]=q+1;
  int*perm=malloc(n*sizeof(int)),*inv=malloc(n*sizeof(int)); long nz; double fl;
  for(int i=0;i<n;i++){perm[i]=i+1;inv[i]=i+1;} stats(n,ptr,row,perm,inv,&nz,&fl);
  printf("natural   nnz(L)=%11ld flops~%.3e\n",nz,fl);
  int64_t ab,ob; fsparse_native_amd_bounded_size(n,q,ptr,row,0,&ab,&ob); int32_t*ar=malloc(ab);
  fsparse_native_amd_bounded(n,q,ptr,row,0,ar,ab,perm,inv); stats(n,ptr,row,perm,inv,&nz,&fl);
  printf("own AMD   nnz(L)=%11ld flops~%.3e\n",nz,fl);
  int*ap=malloc((n+1)*sizeof(int)),*ai=malloc(q*sizeof(int)),*p=malloc(n*sizeof(int));
  for(int i=0;i<=n;i++)ap[i]=ptr[i]-1; for(int i=0;i<q;i++)ai[i]=row[i]-1;
  amd_order(n,ap,ai,p,NULL,NULL); for(int k=0;k<n;k++){perm[k]=p[k]+1;inv[p[k]]=k+1;} stats(n,ptr,row,perm,inv,&nz,&fl);
  printf("SS AMD    nnz(L)=%11ld flops~%.3e\n",nz,fl);
  /* METIS on the off-diagonal graph */
  idx_t nn=n,*xadj=malloc((n+1)*sizeof(idx_t)),*adj=malloc(q*sizeof(idx_t)),*mp=malloc(n*sizeof(idx_t)),*mi=malloc(n*sizeof(idx_t)); idx_t e=0;
  for(int c=0;c<n;c++){xadj[c]=e; for(int t=ptr[c]-1;t<ptr[c+1]-1;t++) if(row[t]-1!=c) adj[e++]=row[t]-1;} xadj[n]=e;
  METIS_NodeND(&nn,xadj,adj,NULL,NULL,mp,mi);
  for(int k=0;k<n;k++){perm[k]=mp[k]+1;inv[mp[k]]=k+1;} stats(n,ptr,row,perm,inv,&nz,&fl);
  printf("METIS ND  nnz(L)=%11ld flops~%.3e\n",nz,fl);
  return 0;}
