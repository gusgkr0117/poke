#include "b-theta4-two.h"
#include <stdlib.h>
#include <string.h>

#ifndef B_THETA4_TWO_CACHE_VECTORS
#define B_THETA4_TWO_CACHE_VECTORS 8
#endif
static unsigned first_bit(unsigned a){unsigned i=0;while(!(a&1u)){a>>=1;i++;}return i;}
static int weight(unsigned a){int w=0;while(a){w+=(int)(a&1u);a>>=1;}return w;}

/* A bounded cache: recursion copies each result to a depth-local scratch
 * vector before making the next call, so eviction cannot invalidate inputs. */
typedef struct {
  b_theta4_two_ctx *Q;
  int evaluation;
  const fq_struct *x;
  const fq_struct *const *shifted;
  fq_struct *unit_scales,*scratch,*scratch_scales,*cache,*cache_scales;
  unsigned cache_index[B_THETA4_TWO_CACHE_VECTORS];
  uint64_t ages[B_THETA4_TWO_CACHE_VECTORS],clock,calls;
  unsigned used;
} b_two_work;
static void work_clear(b_two_work *w){
  b_theta4_ctx *T=w->Q->source;
  b_theta4_vector_free(w->unit_scales,T->g,T->field);
  b_theta4_vector_free(w->scratch,3*T->g*T->count,T->field);
  b_theta4_vector_free(w->scratch_scales,3*T->g,T->field);
  b_theta4_vector_free(w->cache,B_THETA4_TWO_CACHE_VECTORS*T->count,T->field);
  b_theta4_vector_free(w->cache_scales,B_THETA4_TWO_CACHE_VECTORS,T->field);
}
static int work_init(b_two_work *w,b_theta4_two_ctx *Q,int evaluation){
  memset(w,0,sizeof(*w));w->Q=Q;w->evaluation=evaluation;
  b_theta4_ctx *T=Q->source;
  w->unit_scales=b_theta4_vector_new(T->g,T->field);
  w->scratch=b_theta4_vector_new(3*T->g*T->count,T->field);
  w->scratch_scales=b_theta4_vector_new(3*T->g,T->field);
  w->cache=b_theta4_vector_new(B_THETA4_TWO_CACHE_VECTORS*T->count,T->field);
  w->cache_scales=b_theta4_vector_new(B_THETA4_TWO_CACHE_VECTORS,T->field);
  return w->unit_scales&&w->scratch&&w->scratch_scales&&w->cache&&w->cache_scales;
}
static int shift_get(fq_struct *out,fq_t scale,unsigned s,unsigned depth,b_two_work *w){
  b_theta4_ctx *T=w->Q->source;
  if(!w->evaluation&&w->Q->kernel_lift[s]){
    b_theta4_copy(out,w->Q->kernel_lift[s],T);fq_one(scale,T->field);return 1;
  }
  if(w->evaluation&&!s){b_theta4_copy(out,w->x,T);fq_one(scale,T->field);return 1;}
  if(w->evaluation&&weight(s)==1){unsigned i=first_bit(s);
    b_theta4_copy(out,w->shifted[i],T);fq_set(scale,w->unit_scales+i,T->field);return 1;
  }
  for(unsigned c=0;c<w->used;c++)if(w->cache_index[c]==s){
    w->ages[c]=++w->clock;b_theta4_copy(out,w->cache+c*T->count,T);
    fq_set(scale,w->cache_scales+c,T->field);return 1;
  }
  if(depth>=T->g)return 0;
  unsigned i=first_bit(s),j=first_bit(s^(1u<<i));
  unsigned r=s^(1u<<i)^(1u<<j),ri=r|(1u<<i),rj=r|(1u<<j);
  fq_struct *p=w->scratch+3*depth*T->count,*pi=p+T->count,*pj=pi+T->count;
  fq_struct *sp=w->scratch_scales+3*depth,*si=sp+1,*sj=sp+2;
  if(!shift_get(p,sp,r,depth+1,w)||!shift_get(pi,si,ri,depth+1,w)||
     !shift_get(pj,sj,rj,depth+1,w))return 0;
  if(!b_theta4_three_way_add(out,p,w->Q->kernel_lift[1u<<i],
      w->Q->kernel_lift[1u<<j],pi,w->Q->kernel_lift[(1u<<i)|(1u<<j)],pj,T))return 0;
  w->calls++;
  if(fq_is_zero(sp,T->field))return 0;
  fq_mul(scale,si,sj,T->field);fq_div(scale,scale,sp,T->field);
  unsigned slot;
  if(w->used<B_THETA4_TWO_CACHE_VECTORS)slot=w->used++;
  else{slot=0;for(unsigned c=1;c<w->used;c++)if(w->ages[c]<w->ages[slot])slot=c;}
  w->cache_index[slot]=s;w->ages[slot]=++w->clock;
  b_theta4_copy(w->cache+slot*T->count,out,T);fq_set(w->cache_scales+slot,scale,T->field);
  return 1;
}
void b_theta4_two_clear(b_theta4_two_ctx *Q){
  if(Q->source){
    if(Q->kernel_lift)for(size_t s=0;s<Q->source->characters;s++)
      b_theta4_vector_free(Q->kernel_lift[s],Q->source->count,Q->source->field);
    b_theta4_vector_free(Q->target_null,Q->source->count,Q->source->field);
  }
  free(Q->kernel_lift);memset(Q,0,sizeof(*Q));
}
int b_theta4_two_init(b_theta4_two_ctx *Q,b_theta4_ctx *T,
    const fq_struct *const *basis,const fq_struct *const *pairs,const fq_t zeta){
  b_theta4_two_ctx q={0};b_two_work work;int status=0,have_work=0;
  fq_struct *shift=NULL;fq_t scale;fq_init(scale,T->field);q.source=T;
  q.kernel_lift=calloc(T->characters,sizeof(*q.kernel_lift));
  q.target_null=b_theta4_vector_new(T->count,T->field);
  shift=b_theta4_vector_new(T->count,T->field);
  if(!q.kernel_lift||!q.target_null||!shift)goto fail;
  for(unsigned s=0;s<T->characters;s++)if(weight(s)<=2){
    q.kernel_lift[s]=b_theta4_vector_new(T->count,T->field);
    if(!q.kernel_lift[s]){status=0;goto fail;}
    if(s==0)b_theta4_copy(q.kernel_lift[s],T->null,T);
    else if(weight(s)==1){unsigned i=first_bit(s);
      status=b_theta4_good_basis2(q.kernel_lift[s],basis[i],1u<<(2*i),zeta,T);
      if(status!=1)goto fail;
    }else{unsigned i=first_bit(s),j=first_bit(s^(1u<<i));
      status=b_theta4_good_coset2(q.kernel_lift[s],pairs[i*T->g+j],
        q.kernel_lift[1u<<i],q.kernel_lift[1u<<j],1u<<(2*i),zeta,T);
      if(status!=1)goto fail;
    }
  }
  have_work=1;if(!work_init(&work,&q,0)){status=0;goto fail;}
  /* Normal generation and character projection guarantee that at least one
   * binary seed survives. Its single choice is retained for all evaluations. */
  for(q.trace_seed=0;q.trace_seed<T->characters;q.trace_seed++){
    for(unsigned s=0;s<T->characters;s++){
      if(!shift_get(shift,scale,s,0,&work)){status=0;goto fail;}
      b_theta4_two_trace_block(q.target_null,s,q.trace_seed,shift,scale,T);
    }
    if(b_theta4_nonzero(q.target_null,T))break;
  }
  if(q.trace_seed==T->characters){status=0;goto fail;}
  q.init_three_way_calls=work.calls;work_clear(&work);
  b_theta4_vector_free(shift,T->count,T->field);fq_clear(scale,T->field);*Q=q;return 1;
fail:
  if(have_work)work_clear(&work);b_theta4_two_clear(&q);memset(Q,0,sizeof(*Q));
  b_theta4_vector_free(shift,T->count,T->field);fq_clear(scale,T->field);return status;
}
int b_theta4_two_evaluate(fq_struct *out,const fq_struct *x,
    const fq_struct *const *shifted,const fq_t zeta,b_theta4_two_ctx *Q){
  b_theta4_ctx *T=Q->source;fq_struct *result=NULL,*shift=NULL;b_two_work work;
  fq_t one,scale;int ok=0;
  if(!T||!Q->kernel_lift)return 0;
  fq_init(one,T->field);fq_one(one,T->field);fq_init(scale,T->field);
  result=b_theta4_vector_new(T->count,T->field);shift=b_theta4_vector_new(T->count,T->field);
  if(!work_init(&work,Q,1)||!result||!shift)goto done;
  work.x=x;work.shifted=shifted;
  for(unsigned i=0;i<T->g;i++)
    if(!b_theta4_good_coset2_squared_scale(work.unit_scales+i,shifted[i],
      Q->kernel_lift[1u<<i],one,x,1u<<(2*i),zeta,T))goto done;
  for(unsigned s=0;s<T->characters;s++){
    if(!shift_get(shift,scale,s,0,&work))goto done;
    b_theta4_two_trace_block(result,s,Q->trace_seed,shift,scale,T);
  }
  if(!b_theta4_nonzero(result,T))goto done;
  b_theta4_copy(out,result,T);ok=1;
done:
  Q->last_evaluate_three_way_calls=work.calls;work_clear(&work);
  b_theta4_vector_free(result,T->count,T->field);b_theta4_vector_free(shift,T->count,T->field);
  fq_clear(one,T->field);fq_clear(scale,T->field);return ok;
}
int b_theta4_two_evaluate_point(fq_struct *out,const fq_struct *x,
    const fq_t zeta,b_theta4_two_ctx *Q){
  b_theta4_ctx *T=Q->source;if(!T||!Q->kernel_lift)return 0;
  fq_struct *storage=b_theta4_vector_new(T->g*T->count,T->field);
  const fq_struct *shifted[8];int ok=0;if(!storage)return 0;
  for(unsigned i=0;i<T->g;i++){
    shifted[i]=storage+i*T->count;
    if(!b_theta4_add(storage+i*T->count,x,Q->kernel_lift[1u<<i],T))goto done;
  }
  ok=b_theta4_two_evaluate(out,x,shifted,zeta,Q);
done:b_theta4_vector_free(storage,T->g*T->count,T->field);return ok;
}
