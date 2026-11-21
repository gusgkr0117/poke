#include "root-theta-symplectic.h"
#include <string.h>
#include <stdlib.h>
static int mod4(int x){return (x%4+4)%4;}
void root_sp_identity(root_sp_matrix M,unsigned g){memset(M,0,sizeof(root_sp_matrix));for(unsigned j=0;j<2*g;j++)M[j][j]=1;}
int root_sp_valid(const root_sp_matrix M,unsigned g){
 if(!g||g>8)return 0;
 for(unsigned i=0;i<2*g;i++)for(unsigned j=0;j<2*g;j++){
  int x=0;for(unsigned k=0;k<g;k++)x+=M[k][i]*M[g+k][j]-M[g+k][i]*M[k][j];
  int expected=i<g&&j==i+g?1:j<g&&i==j+g?3:0;if(mod4(x)!=expected)return 0;
 }
 return 1;
}
static void matrix_op(root_sp_matrix M,unsigned g,root_sp_op op){
 int i=op.i,j=op.j,c=op.c;for(unsigned k=0;k<2*g;k++){
  int x=M[i][k],y=M[g+i][k];
  switch(op.type){
   case ROOT_SP_FOURIER:if(c>0){M[i][k]=(unsigned char)y;M[g+i][k]=(unsigned char)mod4(-x);}else{M[i][k]=(unsigned char)mod4(-y);M[g+i][k]=(unsigned char)x;}break;
   case ROOT_SP_SWAP:{int z=M[j][k],w=M[g+j][k];M[i][k]=(unsigned char)z;M[j][k]=(unsigned char)x;M[g+i][k]=(unsigned char)w;M[g+j][k]=(unsigned char)y;}break;
   case ROOT_SP_SCALE:M[i][k]=(unsigned char)mod4(c*x);M[g+i][k]=(unsigned char)mod4(c*y);break;
   case ROOT_SP_ADD:M[i][k]=(unsigned char)mod4(x+c*M[j][k]);M[g+j][k]=(unsigned char)mod4(M[g+j][k]-c*y);break;
   case ROOT_SP_LOWER:M[g+i][k]=(unsigned char)mod4(y+c*M[j][k]);if(i!=j)M[g+j][k]=(unsigned char)mod4(M[g+j][k]+c*x);break;
   case ROOT_SP_UPPER:M[i][k]=(unsigned char)mod4(x+c*M[g+j][k]);if(i!=j)M[j][k]=(unsigned char)mod4(M[j][k]+c*y);break;
  }
 }
}
static int push(root_sp_word *W,root_sp_matrix M,int type,int i,int j,int c){
 if(W->count>=ROOT_SP_MAX_OPS)return 0;root_sp_op op={type,i,j,c};W->ops[W->count++]=op;matrix_op(M,W->g,op);return 1;
}
static root_sp_op inverse_op(root_sp_op x){if(x.type!=ROOT_SP_SWAP&&x.type!=ROOT_SP_SCALE)x.c=-x.c;return x;}
void root_sp_word_matrix(root_sp_matrix M,const root_sp_word *W){root_sp_identity(M,W->g);for(unsigned j=0;j<W->count;j++)matrix_op(M,W->g,W->ops[j]);}
int root_sp_inverse(root_sp_matrix out,const root_sp_matrix M,unsigned g){
 if(!root_sp_valid(M,g))return 0;root_sp_matrix N;memset(N,0,sizeof(N));
 for(unsigned i=0;i<g;i++)for(unsigned j=0;j<g;j++){
  N[i][j]=M[g+j][g+i];N[i][g+j]=(unsigned char)mod4(-M[j][g+i]);
  N[g+i][j]=(unsigned char)mod4(-M[g+j][i]);N[g+i][g+j]=M[j][i];
 }
 memcpy(out,N,sizeof(N));return 1;
}
int root_sp_complete_lagrangian(root_sp_matrix out,const root_sp_matrix K,unsigned g){
 if(!g||g>8)return 0;int aug[8][24]={{0}},pivot[8];root_sp_matrix T,M;memset(T,0,sizeof(T));memset(M,0,sizeof(M));
 for(unsigned i=0;i<g;i++)for(unsigned j=0;j<g;j++){
  int s=0;for(unsigned k=0;k<g;k++)s+=K[k][i]*K[g+k][j]-K[g+k][i]*K[k][j];if(mod4(s))return 0;
 }
 for(unsigned i=0;i<g;i++){for(unsigned j=0;j<g;j++){aug[i][j]=mod4(-K[g+j][i]);aug[i][g+j]=mod4(K[j][i]);}aug[i][2*g+i]=1;}
 unsigned rank=0;for(unsigned col=0;col<2*g&&rank<g;col++){
  unsigned row=rank;while(row<g&&!(aug[row][col]&1))row++;if(row==g)continue;
  if(row!=rank)for(unsigned j=0;j<3*g;j++){int s=aug[row][j];aug[row][j]=aug[rank][j];aug[rank][j]=s;}
  if(aug[rank][col]==3)for(unsigned j=0;j<3*g;j++)aug[rank][j]=mod4(3*aug[rank][j]);
  for(unsigned i=0;i<g;i++)if(i!=rank){int s=aug[i][col];for(unsigned j=0;j<3*g;j++)aug[i][j]=mod4(aug[i][j]-s*aug[rank][j]);}
  pivot[rank++]=(int)col;
 }
 if(rank!=g)return 0;for(unsigned i=0;i<g;i++)for(unsigned j=0;j<g;j++)T[pivot[i]][j]=(unsigned char)aug[i][2*g+j];
 int C[8][8]={{0}};for(unsigned i=0;i<g;i++)for(unsigned j=i+1;j<g;j++){
  int s=0;for(unsigned k=0;k<g;k++)s+=T[k][i]*T[g+k][j]-T[g+k][i]*T[k][j];C[i][j]=mod4(s);
 }
 for(unsigned i=0;i<2*g;i++)for(unsigned j=0;j<g;j++){
  M[i][j]=K[i][j]&3;int s=T[i][j];for(unsigned k=0;k<j;k++)s+=K[i][k]*C[k][j];M[i][g+j]=(unsigned char)mod4(s);
 }
 if(!root_sp_valid(M,g))return 0;memcpy(out,M,sizeof(M));return 1;
}
int root_sp_decompose(root_sp_word *out,const root_sp_matrix input,unsigned g){
 if(!root_sp_valid(input,g))return 0;root_sp_matrix M;memcpy(M,input,sizeof(M));root_sp_word R;memset(&R,0,sizeof(R));R.g=g;
 for(unsigned k=0;k<g;k++){
  unsigned pivot=k;int found=0;
  for(unsigned i=k;i<g;i++)if(M[i][k]&1){pivot=i;found=1;break;}
  if(!found)for(unsigned i=k;i<g;i++)if(M[g+i][k]&1){pivot=i;if(!push(&R,M,ROOT_SP_FOURIER,(int)i,0,1))return 0;found=1;break;}
  if(!found)return 0;if(pivot!=k&&!push(&R,M,ROOT_SP_SWAP,(int)k,(int)pivot,1))return 0;
  if(M[k][k]==3&&!push(&R,M,ROOT_SP_SCALE,(int)k,0,3))return 0;
  for(unsigned i=k+1;i<g;i++)if(M[i][k]&&!push(&R,M,ROOT_SP_ADD,(int)i,(int)k,-M[i][k]))return 0;
  for(unsigned i=k+1;i<g;i++)if(M[g+i][k]&&!push(&R,M,ROOT_SP_LOWER,(int)i,(int)k,-M[g+i][k]))return 0;
  if(M[g+k][k]&&!push(&R,M,ROOT_SP_LOWER,(int)k,(int)k,-M[g+k][k]))return 0;
  /* Column k is now e_k. Normalize its symplectic partner while fixing it. */
  if(M[g+k][g+k]!=1)return 0;
  for(unsigned i=k+1;i<g;i++)if(M[g+i][g+k]&&!push(&R,M,ROOT_SP_ADD,(int)k,(int)i,M[g+i][g+k]))return 0;
  for(unsigned i=k+1;i<g;i++)if(M[i][g+k]&&!push(&R,M,ROOT_SP_UPPER,(int)i,(int)k,-M[i][g+k]))return 0;
  if(M[k][g+k]&&!push(&R,M,ROOT_SP_UPPER,(int)k,(int)k,-M[k][g+k]))return 0;
 }
 for(unsigned i=0;i<2*g;i++)for(unsigned j=0;j<2*g;j++)if(M[i][j]!=(i==j))return 0;
 root_sp_word W;memset(&W,0,sizeof(W));W.g=g;W.count=R.count;for(unsigned j=0;j<R.count;j++)W.ops[j]=inverse_op(R.ops[R.count-1-j]);
 root_sp_word_matrix(M,&W);for(unsigned i=0;i<2*g;i++)for(unsigned j=0;j<2*g;j++)if(M[i][j]!=input[i][j])return 0;*out=W;return 1;
}
void root_sp_label(uint32_t *alpha,uint32_t *beta,const root_sp_matrix M,uint32_t a,uint32_t b,unsigned g){
 unsigned char v[16];uint32_t r=0,s=0;for(unsigned i=0;i<g;i++){v[i]=(a>>(2*i))&3u;v[g+i]=(b>>(2*i))&3u;}
 for(unsigned i=0;i<2*g;i++){unsigned t=0;for(unsigned j=0;j<2*g;j++)t+=M[i][j]*v[j];if(i<g)r|=(t&3u)<<(2*i);else s|=(t&3u)<<(2*(i-g));}*alpha=r;*beta=s;
}
int root_theta4_compatibility_sign(int *sign,const fq_struct *P,uint32_t alpha,const b_theta4_ctx *T){
 fq_struct *inv=b_theta4_vector_new(T->count,T->field),*shift=b_theta4_vector_new(T->count,T->field);int ok=0;if(!inv||!shift)goto done;
 fq_t one,c;fq_init(one,T->field);fq_init(c,T->field);fq_one(one,T->field);b_theta4_inverse(inv,P,T);b_theta4_action(shift,inv,alpha,0,one,T);
 if(b_theta4_scale_factor(c,P,shift,T)){if(fq_is_one(c,T->field)){*sign=1;ok=1;}else{fq_neg(c,c,T->field);if(fq_is_one(c,T->field)){*sign=-1;ok=1;}}}
 fq_clear(c,T->field);fq_clear(one,T->field);
done:b_theta4_vector_free(shift,T->count,T->field);b_theta4_vector_free(inv,T->count,T->field);return ok;
}
int root_theta4_compatibility_apply(fq_struct *out,const fq_struct *in,unsigned mask,unsigned g,const fq_ctx_t field){
 if(!g||g>8||(mask>>g))return 0;size_t n=(size_t)1<<(2*g);for(size_t i=0;i<n;i++){
  unsigned odd=0;for(unsigned j=0;j<g;j++)if((mask>>j)&1u)odd^=(unsigned)((i>>(2*j))&1u);
  if(odd)fq_neg(out+i,in+i,field);else fq_set(out+i,in+i,field);
 }return 1;
}
static void fourier(fq_struct *out,const fq_struct *in,unsigned axis,unsigned g,int inverse,const fq_struct *powers,const fq_ctx_t f){
 size_t n=(size_t)1<<(2*g);unsigned shift=2*axis;size_t mask=(size_t)3<<shift;fq_t v,inv4;fq_init(v,f);fq_init(inv4,f);fq_set_ui(inv4,4,f);fq_inv(inv4,inv4,f);
 for(size_t k=0;k<n;k++){unsigned digit=(k>>shift)&3u;fq_zero(out+k,f);for(unsigned j=0;j<4;j++){
  unsigned exponent=(inverse?8u-2u*digit*j:2u*digit*j)&7u;fq_mul(v,in+((k&~mask)|((size_t)j<<shift)),powers+exponent,f);fq_add(out+k,out+k,v,f);
 }if(inverse)fq_mul(out+k,out+k,inv4,f);}
 fq_clear(inv4,f);fq_clear(v,f);
}
static void lower(fq_struct *out,const fq_struct *in,unsigned g,int i,int j,int c,const fq_struct *powers,const fq_ctx_t f){
 size_t n=(size_t)1<<(2*g);for(size_t k=0;k<n;k++){int x=(int)((k>>(2*i))&3u),y=(int)((k>>(2*j))&3u);int exponent=-(i==j?c*x*x:2*c*x*y);exponent=(exponent%8+8)%8;fq_mul(out+k,in+k,powers+exponent,f);}
}
static int apply_one(fq_struct *a,fq_struct *b,const root_sp_op op,unsigned g,const fq_struct *powers,const fq_ctx_t f){
 size_t n=(size_t)1<<(2*g);int i=op.i,j=op.j,c=op.c;if(i<0||(unsigned)i>=g||j<0||(unsigned)j>=g)return 0;
 if(op.type==ROOT_SP_FOURIER){fourier(b,a,(unsigned)i,g,c<0,powers,f);}
 else if(op.type==ROOT_SP_LOWER){lower(b,a,g,i,j,c,powers,f);}
 else if(op.type==ROOT_SP_UPPER){
  /* Only the one or two axes supporting C require a Fourier change;
   * the other Fourier factors commute with L_C and cancel exactly. */
  unsigned axes[2]={(unsigned)i,(unsigned)j},axis_count=i==j?1:2;
  for(unsigned q=0;q<axis_count;q++){fourier(b,a,axes[q],g,0,powers,f);for(size_t k=0;k<n;k++)fq_swap(a+k,b+k,f);}
  lower(b,a,g,i,j,-c,powers,f);for(size_t k=0;k<n;k++)fq_swap(a+k,b+k,f);
  for(unsigned q=0;q<axis_count;q++){fourier(b,a,axes[q],g,1,powers,f);for(size_t k=0;k<n;k++)fq_swap(a+k,b+k,f);}
  return 1;
 }else if(op.type==ROOT_SP_SWAP||op.type==ROOT_SP_SCALE||op.type==ROOT_SP_ADD){
  for(size_t k=0;k<n;k++){int x=(int)((k>>(2*i))&3u),y=(int)((k>>(2*j))&3u);size_t index=k;
   if(op.type==ROOT_SP_SWAP){index&=~(((size_t)3<<(2*i))|((size_t)3<<(2*j)));index|=((size_t)y<<(2*i))|((size_t)x<<(2*j));}
   else{int v=op.type==ROOT_SP_SCALE?mod4(c*x):mod4(x-c*y);index&=~((size_t)3<<(2*i));index|=(size_t)v<<(2*i);}
   fq_set(b+k,a+index,f);
  }
 }else return 0;
 for(size_t k=0;k<n;k++)fq_swap(a+k,b+k,f);return 1;
}
int root_sp_apply(fq_struct *out,const fq_struct *in,const root_sp_word *W,const fq_t eighth_root,int inverse,const fq_ctx_t field){
 if(!W->g||W->g>8||W->count>ROOT_SP_MAX_OPS)return 0;size_t n=(size_t)1<<(2*W->g);int ok=0;
 fq_struct *a=b_theta4_vector_new(n,field),*b=b_theta4_vector_new(n,field),*powers=b_theta4_vector_new(8,field);if(!a||!b||!powers)goto done;
 fq_one(powers,field);for(int j=1;j<8;j++)fq_mul(powers+j,powers+j-1,eighth_root,field);
 fq_t w;fq_init(w,field);fq_one(w,field);fq_neg(w,w,field);if(!fq_equal(powers+4,w,field)){fq_clear(w,field);goto done;}fq_clear(w,field);
 for(size_t k=0;k<n;k++)fq_set(a+k,in+k,field);
 for(unsigned r=0;r<W->count;r++){root_sp_op op=inverse?inverse_op(W->ops[W->count-1-r]):W->ops[r];if(!apply_one(a,b,op,W->g,powers,field))goto done;}
 for(size_t k=0;k<n;k++)fq_set(out+k,a+k,field);ok=1;
done:b_theta4_vector_free(powers,8,field);b_theta4_vector_free(b,n,field);b_theta4_vector_free(a,n,field);return ok;
}
