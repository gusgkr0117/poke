#include "root-theta-label.h"
#include <string.h>
static uint32_t add(uint32_t a,uint32_t b,unsigned g){uint32_t c=0;for(unsigned k=0;k<g;k++)c|=(((a>>(2*k))+(b>>(2*k)))&3u)<<(2*k);return c;}
static uint32_t sub(uint32_t a,uint32_t b,unsigned g){uint32_t c=0;for(unsigned k=0;k<g;k++)c|=(((a>>(2*k))-(b>>(2*k)))&3u)<<(2*k);return c;}
static unsigned dot(uint32_t a,uint32_t b,unsigned g){unsigned c=0;for(unsigned k=0;k<g;k++)c+=((a>>(2*k))&3u)*((b>>(2*k))&3u);return c&3u;}
void root_theta_label_clear(root_theta_label_ctx *L){if(L->ready){b_theta4_vector_free(L->ratios,L->theta->g,L->theta->field);fq_clear(L->zeta,L->theta->field);}memset(L,0,sizeof(*L));}
int root_theta_label_init(root_theta_label_ctx *L,const b_theta4_ctx *T,const fq_t zeta){
 memset(L,0,sizeof(*L));if(!T->g||T->g>8)return 0;unsigned g=T->g;
 fq_t t;fq_init(t,T->field);fq_sqr(t,zeta,T->field);fq_neg(t,t,T->field);int primitive=fq_is_one(t,T->field);fq_clear(t,T->field);if(!primitive)return 0;
 uint32_t pivot=0;while(pivot<T->count&&fq_is_zero(T->null+pivot,T->field))pivot++;if(pivot==T->count)return 0;
 unsigned binary_basis[8]={0},rank=0;uint32_t indices[8];
 for(uint32_t j=0;j<T->count&&rank<g;j++)if(!fq_is_zero(T->null+j,T->field)){
  uint32_t difference=sub(j,pivot,g);unsigned mask=0;for(unsigned k=0;k<g;k++)mask|=((difference>>(2*k))&1u)<<k;
  for(unsigned k=0;k<g;k++)if((mask>>k)&1u){if(binary_basis[k])mask^=binary_basis[k];else{binary_basis[k]=mask;indices[rank++]=j;break;}}
 }
 if(rank<g)return 0;
 int aug[8][16]={{0}};for(unsigned i=0;i<g;i++){uint32_t d=sub(indices[i],pivot,g);for(unsigned k=0;k<g;k++)aug[i][k]=(d>>(2*k))&3;aug[i][g+i]=1;}
 for(unsigned k=0;k<g;k++){
  unsigned row=k;while(row<g&&!(aug[row][k]&1))row++;if(row==g)return 0;
  if(row!=k)for(unsigned j=0;j<2*g;j++){int c=aug[row][j];aug[row][j]=aug[k][j];aug[k][j]=c;}
  if(aug[k][k]==3)for(unsigned j=0;j<2*g;j++)aug[k][j]=(3*aug[k][j])&3;
  for(unsigned row2=0;row2<g;row2++)if(row2!=k){int c=aug[row2][k];for(unsigned j=0;j<2*g;j++)aug[row2][j]=(aug[row2][j]-c*aug[k][j])&3;}
 }
 L->ratios=b_theta4_vector_new(g,T->field);if(!L->ratios)return 0;L->theta=T;L->pivot=pivot;fq_init(L->zeta,T->field);fq_set(L->zeta,zeta,T->field);L->ready=1;
 fq_t inv;fq_init(inv,T->field);fq_inv(inv,T->null+pivot,T->field);
 for(unsigned i=0;i<g;i++){L->index[i]=indices[i];fq_mul(L->ratios+i,T->null+indices[i],inv,T->field);for(unsigned j=0;j<g;j++)L->inverse[i][j]=(unsigned char)aug[i][g+j];}
 fq_clear(inv,T->field);return 1;
}
int root_theta_label(uint32_t *out_alpha,uint32_t *out_beta,const fq_struct *P,const root_theta_label_ctx *L){
 if(!L->ready)return 0;const b_theta4_ctx *T=L->theta;const fq_ctx_struct *f=T->field;unsigned g=T->g;int ok=0;
 fq_t base,rot,neg,scale,tmp;fq_init(base,f);fq_init(rot,f);fq_init(neg,f);fq_init(scale,f);fq_init(tmp,f);fq_struct *phases=b_theta4_vector_new(4,f);if(!phases)goto done;
 for(uint32_t alpha=0;alpha<T->count;alpha++){
  uint32_t anchor=sub(L->pivot,alpha,g);if(fq_is_zero(P+anchor,f))continue;unsigned values[8];int matches=1;
  for(unsigned k=0;k<g;k++){
   uint32_t index=sub(L->index[k],alpha,g);fq_mul(base,P+anchor,L->ratios+k,f);
   if(fq_equal(P+index,base,f)){values[k]=0;continue;}fq_neg(neg,base,f);if(fq_equal(P+index,neg,f)){values[k]=2;continue;}
   fq_mul(rot,base,L->zeta,f);if(fq_equal(P+index,rot,f)){values[k]=1;continue;}fq_neg(neg,rot,f);if(fq_equal(P+index,neg,f)){values[k]=3;continue;}matches=0;break;
  }
  if(!matches)continue;uint32_t beta=0;for(unsigned i=0;i<g;i++){unsigned c=0;for(unsigned j=0;j<g;j++)c+=L->inverse[i][j]*values[j];beta|=(c&3u)<<(2*i);}
  fq_pow_ui(tmp,L->zeta,dot(beta,anchor,g),f);fq_mul(tmp,tmp,T->null+L->pivot,f);fq_div(scale,P+anchor,tmp,f);fq_set(phases,scale,f);for(unsigned i=1;i<4;i++)fq_mul(phases+i,phases+i-1,L->zeta,f);
  for(uint32_t i=0;i<T->count;i++){fq_mul(tmp,phases+dot(beta,i,g),T->null+add(i,alpha,g),f);if(!fq_equal(tmp,P+i,f)){matches=0;break;}}
  if(matches){*out_alpha=alpha;*out_beta=beta;ok=1;break;}
 }
done:
 if(phases)b_theta4_vector_free(phases,4,f);fq_clear(tmp,f);fq_clear(scale,f);fq_clear(neg,f);fq_clear(rot,f);fq_clear(base,f);return ok;
}
