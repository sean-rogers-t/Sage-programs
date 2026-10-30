/* Small driver for the exact F4 API in msolve v0.4.9.
 * Input: nv ng nt; ng lengths; then nt rows (nv exponents numerator denominator).
 * Characteristic zero uses exact GMP arithmetic, not modular reconstruction.
 * A patched core may stop early; the output then consists only of ideal members.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include "libneogb.h"
int main(int argc,char **argv) {
 if(argc!=4 && argc!=5)return 2;
 FILE *in=fopen(argv[1],"r"); if(!in)return 3;
 int32_t nv,ng,nt; if(fscanf(in,"%d%d%d",&nv,&ng,&nt)!=3)return 4;
 int32_t *lens=malloc(ng*sizeof(int32_t)),*exps=malloc(nt*nv*sizeof(int32_t));
 uint32_t prime=strtoul(argv[3],0,10);
 mpz_t **qq=malloc(2*nt*sizeof(mpz_t*));int32_t *ff=malloc(nt*sizeof(int32_t));
 for(int i=0;i<ng;i++)if(fscanf(in,"%d",lens+i)!=1)return 5;
 mpz_t num,den;mpz_inits(num,den,NULL);
 for(int i=0;i<nt;i++) {
  for(int j=0;j<nv;j++)if(fscanf(in,"%d",exps+i*nv+j)!=1)return 6;
  if(gmp_fscanf(in,"%Zd%Zd",num,den)!=2)return 7;
  if(prime) {
   unsigned long a=mpz_fdiv_ui(num,prime),b=mpz_fdiv_ui(den,prime);
   mpz_set_ui(den,b);mpz_set_ui(num,prime);if(!mpz_invert(den,den,num))return 8;
   ff[i]=(a*mpz_get_ui(den))%prime;
  } else {
   qq[2*i]=malloc(sizeof(mpz_t));qq[2*i+1]=malloc(sizeof(mpz_t));
   mpz_init_set(*qq[2*i],num);mpz_init_set(*qq[2*i+1],den);
  }
 }
 fclose(in);
 int32_t nb,*blens=NULL,*bexps=NULL;void *bcf=NULL;
 f4_julia(malloc,&nb,&blens,&bexps,&bcf,lens,exps,prime?(void*)ff:(void*)qq,
          prime,0,0,nv,ng,17,2,0,0,2,0,0,2);
 int32_t *keep=calloc(nb,sizeof(int32_t)),selected=nb;
 if(argc==5) {
  FILE *filter=fopen(argv[4],"r");if(!filter)return 10;
  selected=0;int index;
  while(fscanf(filter,"%d",&index)==1) {
   if(index<0 || index>=nb)return 11;
   if(!keep[index]){keep[index]=1;selected++;}
  }
  fclose(filter);
 } else for(int i=0;i<nb;i++)keep[i]=1;
 printf("F4_TOTAL_MEMBERS=%d\n",nb);fflush(stdout);
 FILE *out=fopen(argv[2],"w");if(!out)return 9;
 fprintf(out,"%d %d\n",nv,selected);int k=0;
 for(int i=0;i<nb;i++) {
  if(!keep[i]){k+=blens[i];continue;}
  fprintf(out,"%d\n",blens[i]);
  for(int j=0;j<blens[i];j++,k++) {
   if(prime)fprintf(out,"%d",((int32_t*)bcf)[k]);
   else gmp_fprintf(out,"%Zd",((mpz_t*)bcf)[k]);
   for(int v=0;v<nv;v++)fprintf(out," %d",bexps[k*nv+v]);
   fputc('\n',out);
  }
 }
 fclose(out);return 0;
}
