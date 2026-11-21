#include "a-kani-fullchain.h"
#include <stdlib.h>
#include <string.h>
void a_kani_fullchain_clear(a_kani_fullchain *S){
 if(!S||!S->ready)return;size_t n=S->theta.count,g=S->theta.g;
 for(size_t j=0;j<S->group_count;j++)if(S->groups[j].field)b_theta4_vector_free(S->groups[j].points,S->groups[j].count*n,S->groups[j].field->field);
 for(size_t j=0;j<S->prime_count;j++)if(S->primes[j].field)b_theta4_vector_free(S->primes[j].pairing4,g*g,S->primes[j].field->field);
 free(S->groups);free(S->primes);fmpz_clear(S->polarization_multiplier);fq_clear(S->zeta8,S->base->field);fq_clear(S->zeta4,S->base->field);b_theta4_clear(&S->theta);memset(S,0,sizeof(*S));
}
static int set_group(a_kani_group *G,const root_extension *E,const fq_struct *const *points,size_t count,size_t n){
 if(!E||(count&&!points)||count>SIZE_MAX/n)return 0;G->field=E;G->count=count;
 if(count){G->points=b_theta4_vector_new(count*n,E->field);if(!G->points)return 0;}
 for(size_t j=0;j<count;j++){if(!points[j])return 0;for(size_t i=0;i<n;i++)fq_set(G->points+j*n+i,points[j]+i,E->field);}return 1;
}
int a_kani_fullchain_init(a_kani_fullchain *S,unsigned g,const root_extension *base,
 const fq_struct *null,const fq_t zeta4,const fq_t zeta8,
 const a_kani_prime_input *primes,size_t pc,const a_theta_payload_input *payloads,size_t gc){
 if(S&&!S->ready)S->last_error=-1;
 if(!S||S->ready||!base||base->relative_degree!=1||!null||!pc||!primes||(gc&&!payloads)||gc>SIZE_MAX-pc||!g||g>8)return 0;
 if(!b_theta4_init(&S->theta,g,null,base->field))return 0;S->base=base;S->ready=1;fq_init(S->zeta4,base->field);fq_init(S->zeta8,base->field);fq_set(S->zeta4,zeta4,base->field);fq_set(S->zeta8,zeta8,base->field);fmpz_init_set_ui(S->polarization_multiplier,1);
 S->groups=calloc(gc+pc,sizeof(*S->groups));S->primes=calloc(pc,sizeof(*S->primes));if(!S->groups||!S->primes)goto fail;S->group_count=gc+pc;S->payload_count=gc;S->prime_count=pc;size_t n=S->theta.count;
 for(size_t j=0;j<gc;j++)if(!set_group(S->groups+j,payloads[j].field,payloads[j].points,payloads[j].count,n))goto fail;
 for(size_t j=0;j<pc;j++){
  const a_kani_prime_input *P=primes+j;a_kani_prime *Q=S->primes+j;
  if(!P->field||!P->exponent||P->ell<2||!P->generators||(P->ell==2&&j)||(P->ell!=2&&(!(P->ell&1)||!P->complements||!P->pairing4)))goto fail;
  for(size_t k=0;k<j;k++)if(primes[k].ell==P->ell)goto fail;
  Q->ell=P->ell;Q->exponent=P->exponent;Q->field=P->field;
  const fq_struct *frame[16];for(unsigned i=0;i<g;i++){frame[i]=P->generators[i];if(P->ell!=2)frame[g+i]=P->complements[i];}
  if(!set_group(S->groups+gc+j,P->field,frame,P->ell==2?g:2*g,n))goto fail;
  if(P->ell!=2){Q->pairing4=b_theta4_vector_new(g*g,P->field->field);if(!Q->pairing4)goto fail;for(unsigned i=0;i<g*g;i++)fq_set(Q->pairing4+i,P->pairing4+i,P->field->field);}
 }
 S->last_error=0;return 1;
fail:a_kani_fullchain_clear(S);S->last_error=-1;return 0;
}
int a_kani_fullchain_next(a_kani_fullchain *S,a_kani_progress progress,void *context){
 if(S)S->last_error=-1;
 if(!S||!S->ready||S->completed==S->prime_count)return 0;size_t n=S->theta.count;unsigned g=S->theta.g;
 a_kani_prime *P=S->primes+S->completed;a_kani_group *frame=S->groups+S->payload_count+S->completed;
 a_theta_payload_input *input=calloc(S->group_count,sizeof(*input));const fq_struct ***ptrs=calloc(S->group_count,sizeof(*ptrs));int ok=0;
 b_theta_balanced binary={0};a_theta_oddbalanced odd={0};fq_struct *pairing4=NULL;b_theta4_ctx target={0};int target_ready=0;
 if(!input||!ptrs)goto done;
 for(size_t j=0;j<S->group_count;j++){
  a_kani_group *G=S->groups+j;ptrs[j]=G->count?calloc(G->count,sizeof(**ptrs)):NULL;if(G->count&&!ptrs[j])goto done;
  for(size_t k=0;k<G->count;k++)ptrs[j][k]=G->points+k*n;input[j].field=G->field;input[j].count=G->count;input[j].points=ptrs[j];
 }
 const fq_struct *K[8],*L[8];for(unsigned i=0;i<g;i++){K[i]=frame->points+i*n;L[i]=P->ell==2?NULL:frame->points+(g+i)*n;}
 if(P->ell==2){
  if(!b_theta_balanced_init(&binary,g,P->exponent,S->base,P->field,S->theta.null,K,input,S->group_count,S->zeta4,S->zeta8)){if(binary.state.last_error==-2)S->last_error=-2;goto done;}
  for(unsigned s=0;s<P->exponent;s++){if(!b_theta_balanced_step(&binary))goto done;if(progress)progress(context,S->completed,P->ell,s+1,P->exponent);}
  if(!b_theta4_init(&target,g,binary.state.theta.null,S->base->field))goto done;target_ready=1;
 }else{
  pairing4=b_theta4_vector_new(g*g,P->field->field);if(!pairing4)goto done;ulong power=fmpz_fdiv_ui(S->polarization_multiplier,P->ell);if(!power)goto done;
  for(unsigned i=0;i<g*g;i++)fq_pow_ui(pairing4+i,P->pairing4+i,power,P->field->field);
  if(!a_theta_oddbalanced_init(&odd,g,P->ell,P->exponent,S->base,P->field,S->theta.null,K,L,pairing4,input,S->group_count)){if(odd.state.last_error==-2)S->last_error=-2;goto done;}
  for(unsigned s=0;s<P->exponent;s++){if(!a_theta_oddbalanced_step(&odd))goto done;if(progress)progress(context,S->completed,P->ell,s+1,P->exponent);}
  if(!b_theta4_init(&target,g,odd.state.theta.null,S->base->field))goto done;target_ready=1;
 }
 /* All operations have succeeded; transfer vector ownership atomically. */
 b_theta4_clear(&S->theta);S->theta=target;memset(&target,0,sizeof(target));target_ready=0;
 for(size_t j=0;j<S->group_count;j++){
  a_kani_group *G=S->groups+j;b_theta4_vector_free(G->points,G->count*n,G->field->field);
  if(P->ell==2){G->points=binary.state.groups[j].points;binary.state.groups[j].points=NULL;binary.state.groups[j].count=0;}
  else{G->points=odd.state.groups[j].points;odd.state.groups[j].points=NULL;odd.state.groups[j].count=0;}
 }
 for(unsigned s=0;s<P->exponent;s++)fmpz_mul_ui(S->polarization_multiplier,S->polarization_multiplier,P->ell);S->completed++;ok=1;
done:
 if(target_ready)b_theta4_clear(&target);a_theta_oddbalanced_clear(&odd);b_theta_balanced_clear(&binary);b_theta4_vector_free(pairing4,g*g,P->field->field);
 if(ptrs){for(size_t j=0;j<S->group_count;j++)free(ptrs[j]);free(ptrs);}free(input);if(ok)S->last_error=0;return ok;
}
