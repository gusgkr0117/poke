#include "b-complete-keygen.h"
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <limits.h>
#include <string.h>
#ifdef INKE_MANUSCRIPT_PARAMETERS
#include <inke_field_constants.h>
#endif

void b_complete_public_clear(b_complete_public *P)
{
    if (!P || !P->ready) return;
    b_product_context_clear(&P->product);
    root_extension_clear(&P->field);
    memset(P, 0, sizeof *P);
}

int b_complete_public_init(b_complete_public *P, a_ec_rng rng, void *context,
                           unsigned long budget)
{
    if (!P || P->ready || !rng) return INKE_SAMPLE_INVALID;
    if (!root_extension_init_rng(&P->field, QUATALG_PINFTY.p, 1, rng, context))
        return INKE_SAMPLE_INVALID;
    int status = b_product_context_init(&P->product, &P->field, rng, context, budget);
    if (status != INKE_SAMPLE_OK) {
        root_extension_clear(&P->field);
        return status;
    }
    P->ready = 1;
    return INKE_SAMPLE_OK;
}
int b_complete_keygen_init(b_complete_keygen *S,b_complete_public *P,
 unsigned long mixed,unsigned long primes,unsigned long norms,unsigned long units){
 return b_complete_keygen_init_with_fast_trials(S,P,mixed,primes,norms,units,
                                               B_FAST_DEFAULT_TRIALS);
}
int b_complete_keygen_init_with_fast_trials(b_complete_keygen *S,b_complete_public *P,
 unsigned long mixed,unsigned long primes,unsigned long norms,unsigned long units,
 unsigned long fast_trials){
 if(!S||S->ready||!P||!P->ready||!fast_trials)return INKE_SAMPLE_INVALID;
 ibz_init(&S->key.q);
#ifdef INKE_PUBLIC_TWO_EXPONENT
 S->key.a=INKE_PUBLIC_TWO_EXPONENT;
#else
 S->key.a=TORSION_PLUS_EVEN_POWER-2;
#endif
 S->public_context=P;
 (void)mixed;(void)primes;(void)norms;(void)units;
 S->fast_trials=fast_trials;S->ready=1;return INKE_SAMPLE_OK;
}
void b_complete_keygen_clear(b_complete_keygen *S){
 if(!S)return;
 if(S->ready)ibz_finalize(&S->key.q);
 memset(S,0,sizeof *S);
}
int b_complete_keygen_step(inke_sk_t *sk,inke_pk_t *pk,b_complete_keygen *S,
 a_kani_progress progress,void *context){
 if(!S||!S->ready)return INKE_SAMPLE_INVALID;
 b_complete_public *P=S->public_context;
 if(S->key.complete)return INKE_SAMPLE_INVALID;
 /* Select q exactly once. Each unsuccessful work batch retains that q;
  * a verified map is retained across all subsequent finishing retries. */
 if(!S->key.has_q){
  S->key.phase=1;int sample=inke_sample_q(&S->key.q,S->key.a);
  if(sample!=INKE_SAMPLE_OK)return sample;
  S->key.has_q=1;
 }
 if(!inke_q_is_admissible(&S->key.q,S->key.a))return INKE_SAMPLE_INVALID;
 if(!S->key.has_conversion){
  S->key.phase=8;inke_conversion_result conversion;
  if(S->fast_counts.norm_candidates>ULONG_MAX-S->fast_trials)return INKE_SAMPLE_INVALID;
  unsigned long limit=S->fast_counts.norm_candidates+S->fast_trials;
  int fast=b_fast_receiver(&conversion,&S->key.q,S->key.a,limit,
         &S->fast_counts,P->product.rng,P->product.rng_context,&S->key.last_stage);
  if(fast==INKE_SAMPLE_OK){
   S->key.conversion=conversion;S->key.has_conversion=1;
   }else return fast;
 }
 a_kani_progress previous=P->product.progress;
 void *previous_context=P->product.progress_context;
 P->product.progress=progress;P->product.progress_context=context;
 S->key.phase=7;
 int status=inke_keygen_from_conversion_with_fallback(sk,pk,&S->key.q,
        &S->key.conversion,b_product_fallback,&P->product,&S->key.last_stage);
 if(status==INKE_SAMPLE_OK){S->key.complete=1;S->key.phase=0;}
 else S->key.finish_failures++;
 P->product.progress=previous;P->product.progress_context=previous_context;
 return status;
}
/* Bounded auxiliary point searches may ask for fresh randomness. Keep the
 * same graph and context until that search succeeds; other errors are final.
 * Intercept the typed status before the PKE API maps it to bool/stage. */
static int decrypt_product(ec_curve_t *target,ec_basis_t *left_images,
 ec_basis_t *right_images,const ec_curve_t *left_curve,
 const ec_curve_t *right_curve,const ec_basis_t *left_two,
 const ec_basis_t *right_two,unsigned int exponent,const ibz_t *q,
 const ec_basis_t *left_points,const ec_basis_t *right_points,
 void *context,int *stage){
 int status;
 do{
  status=b_product_fallback(target,left_images,right_images,left_curve,
      right_curve,left_two,right_two,exponent,q,left_points,right_points,
      context,stage);
 }while(status==INKE_SAMPLE_RETRY);
 return status;
}
int b_complete_decrypt(unsigned char message[32],const inke_ct_t *ct,
 const inke_sk_t *sk,b_complete_public *P,int *stage){
 if(!P||!P->ready)return 0;
 return inke_decrypt_uniform_with_fallback(message,ct,sk,decrypt_product,&P->product,stage);
}
