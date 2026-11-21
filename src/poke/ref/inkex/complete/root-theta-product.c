#include "root-theta-product.h"
#include <stdlib.h>
#include <string.h>
static unsigned digit(unsigned v,unsigned j){return(v>>(2*j))&3;}
static unsigned pairing(unsigned u,unsigned v){return(digit(u,0)*digit(v,2)+digit(u,1)*digit(v,3)+16-digit(u,2)*digit(v,0)+16-digit(u,3)*digit(v,1))&3;}
static unsigned combine(unsigned u,unsigned v,unsigned a,unsigned b){unsigned w=0;for(unsigned j=0;j<4;j++)w|=((a*digit(u,j)+b*digit(v,j))&3)<<(2*j);return w;}
static int primitive(unsigned u){return(u&0x55)!=0;}
static void plane_key(unsigned char key[16],unsigned u,unsigned v){
 unsigned n=0;for(unsigned a=0;a<4;a++)for(unsigned b=0;b<4;b++)key[n++]=(unsigned char)combine(u,v,a,b);
 for(unsigned i=1;i<16;i++){unsigned char x=key[i];unsigned j=i;while(j&&key[j-1]>x){key[j]=key[j-1];j--;}key[j]=x;}
}
int root_theta_product_frames(root_sp_matrix frames[160]){
 unsigned char seen[320][16];size_t found=0,count=0;
 for(unsigned u=1;u<256;u++)if(primitive(u))for(unsigned v=1;v<256;v++)if(pairing(u,v)==1){
  unsigned char key[16],other[16];plane_key(key,u,v);size_t s=0;while(s<found&&memcmp(key,seen[s],16))s++;if(s<found)continue;
  if(found==320)return 0;memcpy(seen[found++],key,16);
  unsigned w=1,z=1;while(w<256&&(!primitive(w)||pairing(u,w)||pairing(v,w)))w++;if(w==256)return 0;
  while(z<256&&(pairing(u,z)||pairing(v,z)||pairing(w,z)!=1))z++;if(z==256)return 0;
  plane_key(other,w,z);if(memcmp(key,other,16)>0)continue;if(count==160)return 0;
  memset(frames[count],0,sizeof(root_sp_matrix));unsigned columns[4]={u,w,v,z};
  for(unsigned j=0;j<4;j++)for(unsigned i=0;i<4;i++)frames[count][i][j]=(unsigned char)digit(columns[j],i);
  if(!root_sp_valid(frames[count],2))return 0;count++;
 }
 return found==320&&count==160;
}
/* Recover the two projective slices, checking every Segre identity. */
static int slices(fq_struct *first,fq_struct *second,const fq_struct *v,const root_extension *E){
 unsigned pivot=0;while(pivot<16&&fq_is_zero(v+pivot,E->field))pivot++;if(pivot==16)return 0;
 unsigned row=pivot%4,column=pivot/4;fq_t x,y;fq_init(x,E->field);fq_init(y,E->field);int ok=1;
 for(unsigned i=0;i<4&&ok;i++)for(unsigned j=0;j<4;j++){
  fq_mul(x,v+i+4*j,v+pivot,E->field);fq_mul(y,v+i+4*column,v+row+4*j,E->field);
  if(!fq_equal(x,y,E->field)){ok=0;break;}
 }
 if(ok)for(unsigned i=0;i<4;i++){fq_set(first+i,v+i+4*column,E->field);fq_set(second+i,v+row+4*i,E->field);}
 fq_clear(y,E->field);fq_clear(x,E->field);return ok;
}
static void heisenberg(fq_struct *out,const fq_struct *in,unsigned a,unsigned b,const root_extension *E){
 fq_t phase;fq_init(phase,E->field);for(unsigned i=0;i<4;i++){fq_pow_ui(phase,E->i,(b*i)&3,E->field);fq_mul(out+i,in+((i+a)&3),phase,E->field);}fq_clear(phase,E->field);
}
/* Internal elliptic basis and theta-characteristic choices are independent
 * on the two tensor axes. Enumerating them separately avoids their product. */
static int elliptic_chart(root_theta_factor *F,root_sp_word *word,unsigned *alpha,unsigned *beta,
 const fq_struct *null,const fq_t eighth,const root_extension *E){
 fq_struct *changed=b_theta4_vector_new(4,E->field),*shifted=b_theta4_vector_new(4,E->field);if(!changed||!shifted){b_theta4_vector_free(shifted,4,E->field);b_theta4_vector_free(changed,4,E->field);return 0;}
 int ok=0;unsigned matrices=0;
 for(unsigned a=0;a<4&&!ok;a++)for(unsigned b=0;b<4&&!ok;b++)for(unsigned c=0;c<4&&!ok;c++)for(unsigned d=0;d<4&&!ok;d++)if(((a*d+16-b*c)&3)==1){
  root_sp_matrix M={{0}};M[0][0]=a;M[0][1]=b;M[1][0]=c;M[1][1]=d;root_sp_word W;matrices++;
  if(!root_sp_decompose(&W,M,1)||!root_sp_apply(changed,null,&W,eighth,0,E->field))goto done;
  for(unsigned x=0;x<4&&!ok;x++)for(unsigned y=0;y<4&&!ok;y++){
   heisenberg(shifted,changed,x,y,E);root_theta_factor candidate;
   if(root_theta_factor_init(&candidate,shifted,1,E)){*F=candidate;*word=W;*alpha=x;*beta=y;ok=1;}
  }
 }
 (void)matrices;
done:b_theta4_vector_free(shifted,4,E->field);b_theta4_vector_free(changed,4,E->field);return ok;
}
void root_theta_product_clear(root_theta_product *P){
 if(P->ready){root_theta_factor_clear(P->factors+1);root_theta_factor_clear(P->factors);fq_clear(P->eighth,P->E->field);}memset(P,0,sizeof(*P));
}
int root_theta_product_init(root_theta_product *P,const fq_struct *null,const fq_t eighth,const root_extension *E){
 if(P->ready)return 0;root_sp_matrix *frames=malloc(160*sizeof(root_sp_matrix));fq_struct *changed=b_theta4_vector_new(16,E->field),*parts=b_theta4_vector_new(8,E->field);int ok=0;
 if(!frames||!changed||!parts)goto done;fq_t check;fq_init(check,E->field);fq_sqr(check,eighth,E->field);int compatible=fq_equal(check,E->i,E->field);fq_clear(check,E->field);if(!compatible||!root_theta_product_frames(frames))goto done;
 for(size_t j=0;j<160;j++){
  root_sp_matrix inverse;root_sp_word common;if(!root_sp_inverse(inverse,frames[j],2)||!root_sp_decompose(&common,inverse,2)||!root_sp_apply(changed,null,&common,eighth,0,E->field))goto done;
  if(!slices(parts,parts+4,changed,E))continue;
  root_theta_product candidate;memset(&candidate,0,sizeof(candidate));candidate.E=E;candidate.common=common;candidate.decompositions_tried=j+1;
  if(!elliptic_chart(candidate.factors,candidate.local,candidate.alpha,candidate.beta,parts,eighth,E))continue;
  if(!elliptic_chart(candidate.factors+1,candidate.local+1,candidate.alpha+1,candidate.beta+1,parts+4,eighth,E)){root_theta_factor_clear(candidate.factors);continue;}
  fq_init(candidate.eighth,E->field);fq_set(candidate.eighth,eighth,E->field);candidate.ready=1;*P=candidate;ok=1;break;
 }
done:b_theta4_vector_free(parts,8,E->field);b_theta4_vector_free(changed,16,E->field);free(frames);return ok;
}
int root_theta_product_point(a_ec_point out[2],const fq_struct *point,const root_theta_product *P){
 if(!P->ready)return 0;const root_extension *E=P->E;fq_struct *v=b_theta4_vector_new(32,E->field);if(!v)return 0;
 a_ec_point result[2];a_ec_init(result,E);a_ec_init(result+1,E);int ok=0;
 if(!root_sp_apply(v,point,&P->common,P->eighth,0,E->field)||!slices(v+16,v+20,v,E))goto done;
 for(unsigned j=0;j<2;j++){
  if(!root_sp_apply(v+24,v+16+4*j,P->local+j,P->eighth,0,E->field))goto done;
  heisenberg(v+28,v+24,P->alpha[j],P->beta[j],E);
  if(!root_theta_factor_point(result+j,v+28,P->factors+j))goto done;
 }
 a_ec_set(out,result,E);a_ec_set(out+1,result+1,E);ok=1;
done:a_ec_clear(result+1,E);a_ec_clear(result,E);b_theta4_vector_free(v,32,E->field);return ok;
}
