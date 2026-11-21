#include "a-theta-oddfields.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
static void up(fq_struct *out,const fq_struct *in,size_t n,const root_field_embedding *L,const root_extension *B,const root_extension *E){for(size_t i=0;i<n;i++){if(B==E)fq_set(out+i,in+i,E->field);else root_field_up(out+i,in+i,L);}}
static int down(fq_struct *out,const fq_struct *in,size_t n,const root_field_embedding *L,const root_extension *B,const root_extension *E){if(B!=E)return root_field_projective_down(out,in,n,L);for(size_t i=0;i<n;i++)fq_set(out+i,in+i,E->field);return 1;}
static slong gcd_s(slong a,slong b){while(b){slong t=a%b;a=b;b=t;}return a;}
void a_theta_oddfields_clear(a_theta_oddfields *S){
 if(!S||!S->ready)return;size_t n=S->theta.count,g=S->theta.g;
 b_theta4_vector_free(S->generators,g*n,S->kernel_field->field);b_theta4_vector_free(S->complements,g*n,S->kernel_field->field);b_theta4_vector_free(S->pairing4,g*g,S->kernel_field->field);
 for(size_t j=0;j<S->group_count;j++){a_theta_odd_payload_group *G=S->groups+j;if(G->field)b_theta4_vector_free(G->points,G->count*n,G->field->field);root_field_embedding_clear(&G->point_up);root_field_embedding_clear(&G->kernel_up);if(G->common_ready)root_extension_clear(&G->common);}
 free(S->groups);root_field_embedding_clear(&S->kernel_embedding);b_theta4_clear(&S->theta);memset(S,0,sizeof(*S));
}
int a_theta_oddfields_init(a_theta_oddfields *S,unsigned g,ulong ell,unsigned exponent,
 const root_extension *B,const root_extension *E,const fq_struct *null,
 const fq_struct *const *K,const fq_struct *const *L,const fq_struct *pairing4,
 const a_theta_payload_input *payloads,size_t groups){
 if(S&&!S->ready)S->last_error=-1;
 if(!S||S->ready||!B||!E||B->relative_degree!=1||!g||g>8||ell<3||!(ell&1)||!exponent||!null||!K||!L||!pairing4||(groups&&!payloads))return 0;
 if(!b_theta4_init(&S->theta,g,null,B->field))return 0;S->ready=1;S->base=B;S->kernel_field=E;S->ell=ell;S->exponent=S->remaining=exponent;
 if(B!=E&&!root_field_embedding_init(&S->kernel_embedding,B,E)){S->last_error=S->kernel_embedding.last_error;goto fail;}
 size_t n=S->theta.count;S->generators=b_theta4_vector_new(g*n,E->field);S->complements=b_theta4_vector_new(g*n,E->field);S->pairing4=b_theta4_vector_new(g*g,E->field);if(!S->generators||!S->complements||!S->pairing4)goto fail;
 for(unsigned j=0;j<g;j++){if(!K[j]||!L[j])goto fail;for(size_t i=0;i<n;i++){fq_set(S->generators+j*n+i,K[j]+i,E->field);fq_set(S->complements+j*n+i,L[j]+i,E->field);}}
 for(unsigned j=0;j<g*g;j++)fq_set(S->pairing4+j,pairing4+j,E->field);
 /* Multiplication commutes with all subsequent quotient maps. Project the
  * complements once instead of repeating this fixed scalar at every stage. */
 if(exponent>1){
  fq_struct *enull=b_theta4_vector_new(n,E->field);b_theta4_ctx T={0};int good=enull!=NULL,ready=0;
  if(good){up(enull,null,n,&S->kernel_embedding,B,E);good=b_theta4_init(&T,g,enull,E->field);ready=good;}
  for(unsigned j=0;good&&j<g;j++)for(unsigned k=1;good&&k<exponent;k++)good=b_theta4_scalar(S->complements+j*n,S->complements+j*n,ell,&T);
  if(ready)b_theta4_clear(&T);b_theta4_vector_free(enull,n,E->field);if(!good)goto fail;
 }
 if(groups){S->groups=calloc(groups,sizeof(*S->groups));if(!S->groups)goto fail;}S->group_count=groups;
 for(size_t j=0;j<groups;j++){
  a_theta_odd_payload_group *G=S->groups+j;const a_theta_payload_input *P=payloads+j;
  if(!P->field||(P->count&&!P->points)||P->count>SIZE_MAX/n)goto fail;G->field=P->field;G->evaluation_field=E;G->count=P->count;
  if(G->count){G->points=b_theta4_vector_new(G->count*n,G->field->field);if(!G->points)goto fail;}
  for(size_t k=0;k<G->count;k++){if(!P->points[k])goto fail;for(size_t i=0;i<n;i++)fq_set(G->points+k*n+i,P->points[k]+i,G->field->field);}
  if(G->field!=E&&G->field!=B){
   slong a=E->relative_degree,b=G->field->relative_degree,d=gcd_s(a,b);if(a/d>WORD_MAX/b)goto fail;
   /* If one point field already contains the other, reuse it literally. */
   if(b%a==0)G->evaluation_field=G->field;
   else if(a%b!=0){if(!root_extension_init_rng(&G->common,fq_ctx_prime(B->field),(a/d)*b,B->rng,B->rng_context)){S->last_error=G->common.last_error;goto fail;}G->common_ready=1;G->evaluation_field=&G->common;}
   if(G->evaluation_field!=E&&!root_field_embedding_init(&G->kernel_up,E,G->evaluation_field)){S->last_error=G->kernel_up.last_error;goto fail;}
   if(G->evaluation_field!=G->field&&!root_field_embedding_init(&G->point_up,G->field,G->evaluation_field)){S->last_error=G->point_up.last_error;goto fail;}
  }
 }
 S->last_error=0;return 1;
fail:{int error=S->last_error;a_theta_oddfields_clear(S);S->last_error=error;}return 0;
}
int a_theta_oddfields_step(a_theta_oddfields *S){
 if(S)S->last_error=-1;
 if(!S||!S->ready||!S->remaining)return 0;
 const root_extension *B=S->base,*E=S->kernel_field;unsigned g=S->theta.g;size_t n=S->theta.count;
 fq_struct *null=b_theta4_vector_new(n,E->field),*small=b_theta4_vector_new(2*g*n,E->field),*next=b_theta4_vector_new(2*g*n,E->field),*targetnull=b_theta4_vector_new(n,B->field);
 fq_struct **groups=S->group_count?calloc(S->group_count,sizeof(*groups)):NULL;
 b_theta4_ctx source={0},target={0};b_theta4_odd_ctx quotient={0};int sr=0,tr=0,ok=0;
 if(!null||!small||!next||!targetnull||(S->group_count&&!groups))goto done;
 up(null,S->theta.null,n,&S->kernel_embedding,B,E);if(!b_theta4_init(&source,g,null,E->field))goto done;sr=1;
 const fq_struct *K[8],*L[8];
 for(unsigned j=0;j<g;j++){
  K[j]=small+j*n;L[j]=small+(g+j)*n;b_theta4_copy(small+j*n,S->generators+j*n,&source);b_theta4_copy(small+(g+j)*n,S->complements+j*n,&source);
  for(unsigned k=1;k<S->remaining;k++)if(!b_theta4_scalar(small+j*n,small+j*n,S->ell,&source))goto done;
 }
 if(!b_theta4_odd_init(&quotient,&source,S->ell,K,L,S->pairing4)||!down(targetnull,quotient.target_null,n,&S->kernel_embedding,B,E)||!b_theta4_init(&target,g,targetnull,B->field))goto done;tr=1;
 for(unsigned j=0;j<g;j++)if(!b_theta4_odd_evaluate(next+j*n,S->generators+j*n,&quotient)||!b_theta4_odd_evaluate(next+(g+j)*n,S->complements+j*n,&quotient))goto done;
 for(size_t k=0;k<S->group_count;k++){
  a_theta_odd_payload_group *G=S->groups+k;if(!G->count)continue;
  const root_extension *C=G->evaluation_field;a_theta_odd_basechange imported={0};b_theta4_odd_ctx *Q=&quotient;
  fq_struct *x=b_theta4_vector_new(2*n,C->field);groups[k]=b_theta4_vector_new(G->count*n,G->field->field);int good=x&&groups[k];
  if(good&&C!=E){good=a_theta_odd_basechange_init(&imported,&quotient,&G->kernel_up);Q=&imported.quotient;}
  for(size_t j=0;good&&j<G->count;j++){
   up(x,G->points+j*n,n,G->point_up.ready?&G->point_up:&S->kernel_embedding,G->field,C);
   good=b_theta4_odd_evaluate(x+n,x,Q);
   if(good)good=down(groups[k]+j*n,x+n,n,G->point_up.ready?&G->point_up:&S->kernel_embedding,G->field,C);
  }
  a_theta_odd_basechange_clear(&imported);b_theta4_vector_free(x,2*n,C->field);if(!good)goto done;
 }
 b_theta4_clear(&S->theta);S->theta=target;memset(&target,0,sizeof(target));tr=0;
 for(size_t i=0;i<g*n;i++){fq_set(S->generators+i,next+i,E->field);fq_set(S->complements+i,next+g*n+i,E->field);}
 for(size_t j=0;j<S->group_count;j++){a_theta_odd_payload_group *G=S->groups+j;b_theta4_vector_free(G->points,G->count*n,G->field->field);G->points=groups[j];groups[j]=NULL;}
 S->remaining--;S->completed++;ok=1;
done:
 b_theta4_odd_clear(&quotient);if(tr)b_theta4_clear(&target);if(sr)b_theta4_clear(&source);
 if(groups){for(size_t j=0;j<S->group_count;j++)b_theta4_vector_free(groups[j],S->groups[j].count*n,S->groups[j].field->field);free(groups);}
 b_theta4_vector_free(targetnull,n,B->field);b_theta4_vector_free(next,2*g*n,E->field);b_theta4_vector_free(small,2*g*n,E->field);b_theta4_vector_free(null,n,E->field);if(ok)S->last_error=0;return ok;
}
