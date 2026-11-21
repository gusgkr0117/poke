#include "b-theta4-odd.h"
#include <string.h>
#include <flint/fmpz.h>
/* Kernel-only normalization is independent of the evaluation point and the
 * selected complement chart. Cache it once when the conservative coefficient
 * storage estimate is <=64MiB; otherwise use the original complete stream. */
static void prepare_cache(b_theta4_odd_ctx *Q){
 if(Q->cache_status)return;Q->cache_status=-1;
#ifdef B_THETA_ODD_DISABLE_CACHE
 return;
#endif
 b_theta4_ctx *T=Q->source;size_t n=T->count,entries=Q->kernel_size;
 size_t limbs=(size_t)((fmpz_bits(fq_ctx_prime(T->field))+FLINT_BITS-1)/FLINT_BITS)+2;
 size_t degree=(size_t)fq_ctx_degree(T->field),limit=((size_t)64*1024*1024)/sizeof(ulong);
 if(!degree||degree>limit/limbs)return;size_t per=degree*limbs;if(n>limit/per||entries>limit/per/(n+1))return;
 fq_struct *points=b_theta4_vector_new(entries*n,T->field),*scales=b_theta4_vector_new(entries,T->field);int ok=points&&scales;ulong digits[8]={0};
 if(ok)b_theta4_copy(points,T->null,T);
 for(size_t k=0;ok&&k<entries;k++){
  if(!b_theta4_odd_kernel_scale(scales+k,points+k*n,Q->ell,T)){ok=0;break;}
  if(k+1==entries)break;b_theta4_copy(points+(k+1)*n,points+k*n,T);
  for(unsigned i=0;i<T->g;i++){
   if(!b_theta4_add(points+(k+1)*n,points+(k+1)*n,Q->basis+i*n,T)){ok=0;break;}
   if(++digits[i]<Q->ell)break;digits[i]=0;
  }
 }
 if(ok){Q->cached_kernel=points;Q->cached_scales=scales;Q->cache_status=1;}
 else{b_theta4_vector_free(scales,entries,T->field);b_theta4_vector_free(points,entries*n,T->field);}
}
int b_theta4_odd_kernel_scale(fq_t scale_ell,const fq_struct *P,ulong ell,b_theta4_ctx *T){
 if(ell<3||!(ell&1u))return 0;
 fq_struct *a=b_theta4_vector_new(T->count,T->field),*b=b_theta4_vector_new(T->count,T->field);int ok=0;
 if(!a||!b)goto done;
 if(!b_theta4_scalar(a,P,ell/2+1,T)||!b_theta4_scalar(b,P,ell/2,T))goto done;
 b_theta4_inverse(b,b,T);ok=b_theta4_scale_factor(scale_ell,a,b,T);
done:b_theta4_vector_free(a,T->count,T->field);b_theta4_vector_free(b,T->count,T->field);return ok;
}
int b_theta4_odd_coset_scale(fq_t scale_ell,const fq_struct *xP,const fq_struct *P,
    const fq_t kernel_scale_ell,const fq_struct *x,ulong ell,b_theta4_ctx *T){
 if(ell<3||!(ell&1u)||fq_is_zero(kernel_scale_ell,T->field))return 0;
 fq_struct *a=b_theta4_vector_new(T->count,T->field),*b=b_theta4_vector_new(T->count,T->field);int ok=0;fq_t scale,power;
 fq_init(scale,T->field);fq_init(power,T->field);if(!a||!b)goto done;
 if(!b_theta4_scalar_add(a,b,ell,xP,P,x,T)||!b_theta4_scale_factor(scale,a,x,T))goto done;
 fq_pow_ui(power,kernel_scale_ell,ell-1,T->field);fq_div(scale,scale,power,T->field);fq_set(scale_ell,scale,T->field);ok=1;
done:b_theta4_vector_free(a,T->count,T->field);b_theta4_vector_free(b,T->count,T->field);fq_clear(scale,T->field);fq_clear(power,T->field);return ok;
}
int b_theta4_odd_term(fq_struct *out,const fq_struct *xP,const fq_struct *P,
    const fq_struct *x,ulong ell,b_theta4_ctx *T){
 fq_t k,c;fq_init(k,T->field);fq_init(c,T->field);int ok=0;
 if(!b_theta4_odd_kernel_scale(k,P,ell,T)||!b_theta4_odd_coset_scale(c,xP,P,k,x,ell,T))goto done;
 for(size_t j=0;j<T->count;j++){fq_pow_ui(out+j,xP+j,ell,T->field);fq_mul(out+j,out+j,c,T->field);}ok=1;
done:fq_clear(k,T->field);fq_clear(c,T->field);return ok;
}

/* Returns 1 even for the zero trace: chart selection must distinguish that
 * legitimate exceptional section from an arithmetic/input failure. */
static int odd_trace(fq_struct *out,const fq_struct *x,b_theta4_odd_ctx *Q){
 b_theta4_ctx *T=Q->source;size_t n=T->count;
 prepare_cache(Q);
 fq_struct *v=b_theta4_vector_new(5*n,T->field);int ok=0;ulong digits[8]={0};
 fq_t coefficient,coset;fq_init(coefficient,T->field);fq_init(coset,T->field);fq_one(coefficient,T->field);
 if(!v)goto done;
 fq_struct *base=v,*P=v+n,*xP=v+2*n,*term=v+3*n,*sum=v+4*n;
 if(Q->selected_chart){if(!b_theta4_add(base,x,Q->translate,T))goto done;}
 else b_theta4_copy(base,x,T);
 b_theta4_copy(P,T->null,T);b_theta4_copy(xP,base,T);
 for(size_t j=0;j<n;j++)fq_zero(sum+j,T->field);
 for(size_t k=0;k<Q->kernel_size;k++){
  if(Q->cache_status==1){
   if(!b_theta4_odd_coset_scale(coset,xP,Q->cached_kernel+k*n,Q->cached_scales+k,base,Q->ell,T))goto done;
   for(size_t j=0;j<n;j++){fq_pow_ui(term+j,xP+j,Q->ell,T->field);fq_mul(term+j,term+j,coset,T->field);}
  }else if(!b_theta4_odd_term(term,xP,P,base,Q->ell,T))goto done;
  for(size_t j=0;j<n;j++){
   fq_mul(term+j,term+j,coefficient,T->field);fq_add(sum+j,sum+j,term+j,T->field);
  }
  if(k+1==Q->kernel_size)break;
  for(unsigned i=0;i<T->g;i++){
   if((Q->cache_status!=1&&!b_theta4_add(P,P,Q->basis+i*n,T))||!b_theta4_add(xP,xP,Q->basis+i*n,T))goto done;
   fq_mul(coefficient,coefficient,Q->weights+i,T->field);
   if(++digits[i]<Q->ell)break;
   digits[i]=0;
  }
 }
 b_theta4_copy(out,sum,T);ok=1;
done:b_theta4_vector_free(v,5*n,T->field);fq_clear(coset,T->field);fq_clear(coefficient,T->field);return ok;
}
void b_theta4_odd_clear(b_theta4_odd_ctx *Q){
 if(Q->source){b_theta4_ctx *T=Q->source;
  b_theta4_vector_free(Q->basis,T->g*T->count,T->field);
  b_theta4_vector_free(Q->translate,T->count,T->field);
  b_theta4_vector_free(Q->weights,T->g,T->field);
  b_theta4_vector_free(Q->target_null,T->count,T->field);
  b_theta4_vector_free(Q->cached_scales,Q->kernel_size,T->field);
  if(Q->cached_kernel)b_theta4_vector_free(Q->cached_kernel,Q->kernel_size*T->count,T->field);
 }
 memset(Q,0,sizeof(*Q));
}
int b_theta4_odd_init(b_theta4_odd_ctx *Q,b_theta4_ctx *T,ulong ell,
 const fq_struct *const *kernel,const fq_struct *const *complement,const fq_struct *pairing4){
 b_theta4_odd_ctx q={0};ulong digits[8]={0};int ok=0;
 if(ell<3||!(ell&1u)){memset(Q,0,sizeof(*Q));return 0;}
 q.source=T;q.ell=ell;q.kernel_size=1;
 for(unsigned i=0;i<T->g;i++){
  if(q.kernel_size>SIZE_MAX/ell)goto done;q.kernel_size*=ell;
 }
 q.basis=b_theta4_vector_new(T->g*T->count,T->field);q.translate=b_theta4_vector_new(T->count,T->field);
 q.weights=b_theta4_vector_new(T->g,T->field);q.target_null=b_theta4_vector_new(T->count,T->field);
 if(!q.basis||!q.translate||!q.weights||!q.target_null)goto done;
 b_theta4_copy(q.translate,T->null,T);
 for(unsigned i=0;i<T->g;i++){b_theta4_copy(q.basis+i*T->count,kernel[i],T);fq_one(q.weights+i,T->field);}
 for(q.selected_chart=0;q.selected_chart<q.kernel_size;q.selected_chart++){
  if(!odd_trace(q.target_null,T->null,&q))goto done;
  if(b_theta4_nonzero(q.target_null,T)){*Q=q;ok=1;break;}
  if(q.selected_chart+1==q.kernel_size)break;
  for(unsigned j=0;j<T->g;j++){
   if(!b_theta4_add(q.translate,q.translate,complement[j],T))goto done;
   for(unsigned i=0;i<T->g;i++)fq_mul(q.weights+i,q.weights+i,pairing4+i*T->g+j,T->field);
   if(++digits[j]<ell)break;digits[j]=0;
  }
 }
done:if(!ok){b_theta4_odd_clear(&q);memset(Q,0,sizeof(*Q));}return ok;
}
int b_theta4_odd_evaluate(fq_struct *out,const fq_struct *x,b_theta4_odd_ctx *Q){
 if(!Q->source||!Q->basis)return 0;
 b_theta4_ctx *T=Q->source;fq_struct *v=b_theta4_vector_new(T->count,T->field);int ok=0;
 if(v&&odd_trace(v,x,Q)&&b_theta4_nonzero(v,T)){b_theta4_copy(out,v,T);ok=1;}
 b_theta4_vector_free(v,T->count,T->field);return ok;
}
