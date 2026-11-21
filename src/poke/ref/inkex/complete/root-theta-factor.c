#include "root-theta-factor.h"
#include <string.h>
static int slice(fq_struct *out,const fq_struct *P,size_t count,const root_extension *E){
 size_t pivot=0;while(pivot<count&&fq_is_zero(P+pivot,E->field))pivot++;if(pivot==count)return 0;
 size_t column=pivot/4;unsigned row=(unsigned)(pivot%4);fq_t a,b;fq_init(a,E->field);fq_init(b,E->field);int ok=1;
 for(size_t j=0;j<count/4&&ok;j++)for(unsigned i=0;i<4;i++){
  fq_mul(a,P+4*j+i,P+pivot,E->field);fq_mul(b,P+4*column+i,P+4*j+row,E->field);if(!fq_equal(a,b,E->field)){ok=0;break;}
 }
 if(ok)for(unsigned i=0;i<4;i++)fq_set(out+i,P+4*column+i,E->field);fq_clear(b,E->field);fq_clear(a,E->field);return ok;
}
int root_theta_factor_init(root_theta_factor *F,const fq_struct *O,unsigned g,const root_extension *E){
 memset(F,0,sizeof(*F));if(!g||g>8)return 0;fq_struct *v=b_theta4_vector_new(4,E->field);if(!v)return 0;size_t count=(size_t)1<<(2*g);
 int ok=slice(v,O,count,E);if(ok){ok=root_ec_theta4_from_null(&F->elliptic,v,E);if(!ok)root_ec_theta4_clear(&F->elliptic);}b_theta4_vector_free(v,4,E->field);
 if(ok){F->E=E;F->g=g;F->count=count;F->ready=1;}return ok;
}
void root_theta_factor_clear(root_theta_factor *F){if(F->ready)root_ec_theta4_clear(&F->elliptic);memset(F,0,sizeof(*F));}
int root_theta_factor_point(a_ec_point *out,const fq_struct *P,const root_theta_factor *F){
 if(!F->ready)return 0;fq_struct *v=b_theta4_vector_new(4,F->E->field);if(!v)return 0;
 int ok=slice(v,P,F->count,F->E)&&root_ec_theta4_unmap(out,v,&F->elliptic);b_theta4_vector_free(v,4,F->E->field);return ok;
}
