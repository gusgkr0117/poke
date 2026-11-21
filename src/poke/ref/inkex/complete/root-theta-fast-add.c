/* An experimental replacement translation unit for b-theta4.c.  Keep
 * its complete addition as a reference/fallback while validating the
 * diagonal Walsh fastpath independently before merging it upstream. */
#include "root-extension.h"
/* The production replacement uses explicit uniform randomness for its
 * generic-context square roots, never the library's fixed-seed search. */
#define fq_sqrt root_fq_sqrt
#define b_theta4_add b_theta4_add_reference
#define B_THETA4_DIFF_IMPL root_theta_diff_reference
#define B_THETA4_THREE_IMPL root_theta_three_reference
#include "b-theta4.c"
#undef b_theta4_add
#undef fq_sqrt

unsigned long root_theta_fast_add_full_support=0;
unsigned long root_theta_fast_add_one_column=0;
unsigned long root_theta_fast_add_fallback=0;
unsigned long root_theta_fast_add_sparse=0;
unsigned long root_theta_fast_add_residual_rows=0;
unsigned long root_theta_fast_add_residual_blocks=0;
unsigned long root_theta_fast_add_residual_checks=0;
unsigned long root_theta_fast_known_diagonal_rows=0,root_theta_fast_known_residual_rows=0;
unsigned long root_theta_fast_known_residual_blocks=0;
unsigned long root_theta_fast_known_affine_checks=0;

/* The inverse Riemann transform is M_ij=2^g U_i V_j. Here V is known
 * with its affine scale, so each nonzero V_j recovers U_i with that SAME
 * scale. Probe whole diagonals before the complete fixed-column fallback.
 * Every accepted coordinate is divided by 2^g V_j; no projective
 * normalization is introduced into good-lift computations. */
static int known_factor(fq_struct *out,const fq_struct *V,
 const fq_struct *P1,const fq_struct *P2,const fq_struct *Q1,const fq_struct *Q2,
 const fq_struct *denA,const fq_struct *denB,b_theta4_chart *charts,b_theta4_ctx *T){
 fq_struct *result=b_theta4_vector_new(T->count,T->field),*block=b_theta4_vector_new(T->characters,T->field),*scratch=b_theta4_vector_new(2*T->characters,T->field);
 unsigned char *done=calloc(T->count,1);fq_t den,sum;fq_init(den,T->field);fq_init(sum,T->field);int ok=0;size_t pivot=0,remaining=T->count;while(pivot<T->count&&fq_is_zero(V+pivot,T->field))pivot++;
 if(!result||!block||!scratch||!done||pivot==T->count)goto finish;
#ifndef ROOT_FAST_FORCE_RESIDUAL
 for(unsigned pass=0;pass<=2*T->g&&remaining;pass++){
  uint32_t delta=pass?((pass<=T->g?1u:3u)<<(2*((pass-1)%T->g))):0;
  for(unsigned q=0;q<T->characters;q++){
   uint32_t base=emb(q,T->g);int need=0;for(unsigned s=0;s<T->characters;s++){uint32_t i=base+2*emb(s,T->g),j=addidx(i,delta,T->g);if(!done[i]&&!fq_is_zero(V+j,T->field)){need=1;break;}}if(!need)continue;
   if(!riemann_all(block,base,addidx(base,delta,T->g),P1,P2,Q1,Q2,denA,denB,charts,scratch,scratch+T->characters,T))goto finish;walsh(block,T->characters,T->field);
   for(unsigned s=0;s<T->characters;s++){uint32_t i=base+2*emb(s,T->g),j=addidx(i,delta,T->g);if(done[i]||fq_is_zero(V+j,T->field))continue;fq_mul_ui(den,V+j,(ulong)T->characters,T->field);fq_inv(den,den,T->field);fq_mul(result+i,block+s,den,T->field);done[i]=1;remaining--;root_theta_fast_known_diagonal_rows++;}
  }
 }
#endif
 /* Cover an entire parity class per inverse Walsh transform when possible.
  * For the selected unresolved row, delta sends its partner to the known
  * nonzero pivot. Thus every block makes progress without a support promise. */
 for(unsigned q=0;q<T->characters&&remaining;q++){
  uint32_t base=emb(q,T->g);
  for(;;){
   unsigned first=0;while(first<T->characters&&done[base+2*emb(first,T->g)])first++;if(first==T->characters)break;
   uint32_t selected=base+2*emb(first,T->g),delta=addidx((uint32_t)pivot,negidx(selected,T->g),T->g);
   if(!riemann_all(block,base,addidx(base,delta,T->g),P1,P2,Q1,Q2,denA,denB,charts,scratch,scratch+T->characters,T))goto finish;walsh(block,T->characters,T->field);root_theta_fast_known_residual_blocks++;
   for(unsigned s=0;s<T->characters;s++){
    uint32_t i=base+2*emb(s,T->g),j=addidx(i,delta,T->g);if(done[i]||fq_is_zero(V+j,T->field))continue;
    fq_mul_ui(den,V+j,(ulong)T->characters,T->field);fq_inv(den,den,T->field);fq_mul(result+i,block+s,den,T->field);done[i]=1;remaining--;root_theta_fast_known_residual_rows++;
   }
   if(!done[selected])goto finish;
  }
 }
 if(!b_theta4_nonzero(result,T))goto finish;b_theta4_copy(out,result,T);ok=1;
finish:free(done);fq_clear(sum,T->field);fq_clear(den,T->field);b_theta4_vector_free(scratch,2*T->characters,T->field);b_theta4_vector_free(block,T->characters,T->field);b_theta4_vector_free(result,T->count,T->field);return ok;
}
int b_theta4_diff_add(fq_struct *out,const fq_struct *P,const fq_struct *Q,const fq_struct *PmQ,b_theta4_ctx *T){
#ifdef ROOT_FAST_VERIFY_KNOWN
 fq_struct *reference=NULL;int expected=0;if(T->g<=3){reference=b_theta4_vector_new(T->count,T->field);if(!reference)return 0;expected=root_theta_diff_reference(reference,P,Q,PmQ,T);}
#endif
 int ok=known_factor(out,PmQ,P,P,Q,Q,T->null,T->null,T->charts,T);
#ifdef ROOT_FAST_VERIFY_KNOWN
 if(reference){if(ok!=expected)abort();if(ok)for(size_t i=0;i<T->count;i++)if(!fq_equal(out+i,reference+i,T->field))abort();root_theta_fast_known_affine_checks++;b_theta4_vector_free(reference,T->count,T->field);}
#endif
 return ok;
}
int b_theta4_three_way_add(fq_struct *out,const fq_struct *P,const fq_struct *Q,const fq_struct *R,const fq_struct *PQ,const fq_struct *QR,const fq_struct *PR,b_theta4_ctx *T){
#ifdef ROOT_FAST_VERIFY_KNOWN
 fq_struct *reference=NULL;int expected=0;if(T->g<=3){reference=b_theta4_vector_new(T->count,T->field);if(!reference)return 0;expected=root_theta_three_reference(reference,P,Q,R,PQ,QR,PR,T);}
#endif
 fq_struct *n=b_theta4_vector_new(3*T->count,T->field);b_theta4_chart *charts=calloc(T->count,sizeof(*charts));if(!n||!charts){b_theta4_vector_free(n,3*T->count,T->field);free(charts);return 0;}for(size_t i=0;i<T->count;i++)fq_init(charts[i].inverse_denominator,T->field);b_theta4_inverse(n,P,T);b_theta4_inverse(n+T->count,Q,T);b_theta4_inverse(n+2*T->count,R,T);
 int ok=known_factor(out,n,T->null,QR,PR,PQ,n+T->count,n+2*T->count,charts,T);for(size_t i=0;i<T->count;i++)fq_clear(charts[i].inverse_denominator,T->field);free(charts);b_theta4_vector_free(n,3*T->count,T->field);
#ifdef ROOT_FAST_VERIFY_KNOWN
 if(reference){if(ok!=expected)abort();if(ok)for(size_t i=0;i<T->count;i++)if(!fq_equal(out+i,reference+i,T->field))abort();root_theta_fast_known_affine_checks++;b_theta4_vector_free(reference,T->count,T->field);}
#endif
 return ok;
}

/* out_i=2^g U_i V_(i+delta), where U and V are ONE consistent pair
 * of affine lifts of P+Q and P-Q. All offsets have the same scalar. */
static int outer_diagonal_prefix(fq_struct *out,uint32_t delta,unsigned parity_count,
    const fq_struct *P,const fq_struct *Q,fq_struct *block,
    fq_struct *scratch,b_theta4_ctx *T){
 for(unsigned q=0;q<parity_count;q++){
  uint32_t base=emb(q,T->g);
  if(!riemann_all(block,base,addidx(base,delta,T->g),P,P,Q,Q,
       T->null,T->null,T->charts,scratch,scratch+T->characters,T))return 0;
  walsh(block,T->characters,T->field);
  for(unsigned s=0;s<T->characters;s++)fq_set(out+base+2*emb(s,T->g),block+s,T->field);
 }
 return 1;
}
static int outer_diagonal(fq_struct *out,uint32_t delta,
    const fq_struct *P,const fq_struct *Q,fq_struct *block,
    fq_struct *scratch,b_theta4_ctx *T){
 return outer_diagonal_prefix(out,delta,T->characters,P,Q,block,scratch,T);
}

/* Recover one consistently scaled rank-one factor on whatever support
 * is reached by sparse diagonals. Missing rows use the certified pivot
 * column; no connectivity or generic-support hypothesis is assumed. */
static int sparse_factor(fq_struct *answer,const fq_struct *D,
    const fq_struct *P,const fq_struct *Q,fq_struct *block,
    fq_struct *scratch,b_theta4_ctx *T){
 size_t n=T->count,pivot=n,first=n;
 fq_struct *edges=b_theta4_vector_new(T->g*n,T->field),*columns=b_theta4_vector_new(n,T->field);
 unsigned char *known=calloc(2*n,1);size_t *queue=malloc(2*n*sizeof(*queue));
 fq_t inverse,value;fq_init(inverse,T->field);fq_init(value,T->field);int ok=0;
 if(!edges||!columns||!known||!queue)goto done;
 for(size_t i=0;i<n;i++)if(!fq_is_zero(D+i,T->field)){first=pivot=i;break;}
 for(unsigned axis=0;axis<T->g;axis++){
  uint32_t delta=1u<<(2*axis);
  if(!outer_diagonal(edges+axis*n,delta,P,Q,block,scratch,T))goto done;
  if(pivot==n)for(size_t i=0;i<n;i++)if(!fq_is_zero(edges+axis*n+i,T->field)){
   first=i;pivot=addidx((uint32_t)i,delta,T->g);fq_set(answer+first,edges+axis*n+i,T->field);break;
  }
 }
 if(pivot==n){ok=-1;goto done;}
 if(first==pivot)fq_set(answer+first,D+first,T->field);
 fq_one(columns+pivot,T->field);known[first]=known[n+pivot]=1;
 size_t head=0,tail=0;queue[tail++]=first;queue[tail++]=n+pivot;
 while(head<tail){
  size_t node=queue[head++],i=node<n?node:node-n;int row=node<n;
  fq_inv(inverse,row?answer+i:columns+i,T->field);
  for(unsigned axis=0;axis<=T->g;axis++){
   uint32_t delta=axis?(1u<<(2*(axis-1))):0;
   size_t other=row?addidx((uint32_t)i,delta,T->g):addidx((uint32_t)i,negidx(delta,T->g),T->g);
   size_t neighbor=row?n+other:other;if(known[neighbor])continue;
   const fq_struct *entry=axis?edges+(axis-1)*n+(row?i:other):D+i;
   fq_mul(value,entry,inverse,T->field);
   fq_set(row?columns+other:answer+other,value,T->field);known[neighbor]=1;
   if(!fq_is_zero(value,T->field))queue[tail++]=neighbor;
  }
 }
 unsigned long residual=0;
 for(unsigned q=0;q<T->characters;q++){
  uint32_t base=emb(q,T->g);
  for(;;){
   unsigned first_unknown=0;while(first_unknown<T->characters&&known[base+2*emb(first_unknown,T->g)])first_unknown++;if(first_unknown==T->characters)break;
   uint32_t selected=base+2*emb(first_unknown,T->g),delta=addidx((uint32_t)pivot,negidx(selected,T->g),T->g);
   if(!riemann_all(block,base,addidx(base,delta,T->g),P,P,Q,Q,T->null,T->null,T->charts,scratch,scratch+T->characters,T))goto done;
   walsh(block,T->characters,T->field);root_theta_fast_add_residual_blocks++;
   for(unsigned s=0;s<T->characters;s++){
    uint32_t i=base+2*emb(s,T->g),j=addidx(i,delta,T->g);
    if(known[i]||!known[n+j]||fq_is_zero(columns+j,T->field))continue;
    fq_div(answer+i,block+s,columns+j,T->field);known[i]=1;residual++;
#ifdef ROOT_FAST_VERIFY_KNOWN
    if(T->g<=3){if(!riemann_sum(value,i,(uint32_t)pivot,P,Q,scratch,scratch+T->characters,T)||!fq_equal(value,answer+i,T->field))abort();root_theta_fast_add_residual_checks++;}
#endif
   }
   if(!known[selected])goto done;
  }
 }
 if(!b_theta4_nonzero(answer,T))goto done;
 root_theta_fast_add_sparse++;root_theta_fast_add_residual_rows+=residual;
 if(residual)root_theta_fast_add_one_column++;ok=1;
done:fq_clear(value,T->field);fq_clear(inverse,T->field);free(queue);free(known);b_theta4_vector_free(columns,n,T->field);b_theta4_vector_free(edges,T->g*n,T->field);return ok;
}

int b_theta4_add(fq_struct *out,const fq_struct *P,const fq_struct *Q,b_theta4_ctx *T){
 fq_struct *D=b_theta4_vector_new(T->count,T->field),*edge=b_theta4_vector_new(T->count,T->field),*answer=b_theta4_vector_new(T->count,T->field),*block=b_theta4_vector_new(T->characters,T->field),*scratch=b_theta4_vector_new(2*T->characters,T->field);
 fq_t ratio;fq_init(ratio,T->field);int ok=0;size_t nonzero=0;
 if(!D||!edge||!answer||!block||!scratch)goto done;
 if(!outer_diagonal(D,0,P,Q,block,scratch,T))goto done;
 for(size_t i=0;i<T->count;i++)if(!fq_is_zero(D+i,T->field))nonzero++;
 if(nonzero==T->count){
  /* A spanning tree adding one positive unit per tensor axis.  Since
   * every diagonal is nonzero, both U and V have full support and
   * A_child=A_parent D_child/M_(parent,child)=U_child/U_0. */
  fq_one(answer,T->field);size_t stride=1;
  for(unsigned axis=0;axis<T->g;axis++,stride*=4){
   /* This spanning tree reads only parents whose higher digits are zero.
    * Their parity representatives therefore use the first axis+1 bits. */
   if(!outer_diagonal_prefix(edge,1u<<(2*axis),1u<<(axis+1),P,Q,block,scratch,T))goto done;
   for(unsigned digit=1;digit<4;digit++)for(size_t low=0;low<stride;low++){
    size_t child=digit*stride+low,parent=child-stride;
    if(fq_is_zero(edge+parent,T->field))goto done;
    fq_div(ratio,D+child,edge+parent,T->field);fq_mul(answer+child,answer+parent,ratio,T->field);
   }
  }
  root_theta_fast_add_full_support++;
 }else{
  int sparse=sparse_factor(answer,D,P,Q,block,scratch,T);
  if(!sparse)goto done;
  if(sparse<0){root_theta_fast_add_fallback++;if(!b_theta4_add_reference(answer,P,Q,T))goto done;}
 }
 b_theta4_copy(out,answer,T);ok=1;
done:fq_clear(ratio,T->field);b_theta4_vector_free(scratch,2*T->characters,T->field);b_theta4_vector_free(block,T->characters,T->field);b_theta4_vector_free(answer,T->count,T->field);b_theta4_vector_free(edge,T->count,T->field);b_theta4_vector_free(D,T->count,T->field);return ok;
}
