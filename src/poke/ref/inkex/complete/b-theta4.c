#include "b-theta4.h"
#include <stdlib.h>
#include <string.h>

/* Formula source: Dequay--Lubicz, arXiv:2504.21058v2, Theorem 3 Eq.18.
 * Index derivation is recorded in dimension8-geometry-plan.md. This is a
 * direct C implementation of the identity, not a binding to Sage/Python. */
static unsigned parity(uint32_t x) {
  unsigned r = 0;
  while (x) {
    r ^= x & 1u;
    x >>= 1;
  }
  return r;
}
static uint32_t addidx(uint32_t a, uint32_t b, unsigned g) {
  uint32_t r = 0;
  unsigned j;
  for (j = 0; j < g; j++)
    r |= (((a >> (2 * j)) + (b >> (2 * j))) & 3u) << (2 * j);
  return r;
}
static uint32_t negidx(uint32_t a, unsigned g) {
  uint32_t r = 0;
  unsigned j;
  for (j = 0; j < g; j++)
    r |= ((4u - ((a >> (2 * j)) & 3u)) & 3u) << (2 * j);
  return r;
}
static uint32_t emb(uint32_t a, unsigned g) {
  uint32_t r = 0;
  unsigned j;
  for (j = 0; j < g; j++)
    r |= ((a >> j) & 1u) << (2 * j);
  return r;
}
static uint32_t lowbits(uint32_t a, unsigned g) {
  uint32_t r = 0;
  unsigned j;
  for (j = 0; j < g; j++)
    r |= ((a >> (2 * j)) & 1u) << j;
  return r;
}
static uint32_t halfidx(uint32_t i, uint32_t j, uint32_t k, uint32_t l,
                        unsigned g) {
  uint32_t r = 0;
  unsigned s;
  for (s = 0; s < g; s++) {
    int v = -(int)((i >> (2 * s)) & 3u) + (int)((j >> (2 * s)) & 3u) +
            (int)((k >> (2 * s)) & 3u) + (int)((l >> (2 * s)) & 3u);
    r |= (uint32_t)(((v / 2) % 4 + 4) % 4) << (2 * s);
  }
  return r;
}
fq_struct *b_theta4_vector_new(size_t n, const fq_ctx_t field) {
  size_t j;
  fq_struct *v = (fq_struct *)calloc(n, sizeof(*v));
  if (!v)
    return NULL;
  for (j = 0; j < n; j++)
    fq_init(v + j, field);
  return v;
}
void b_theta4_vector_free(fq_struct *v, size_t n, const fq_ctx_t field) {
  size_t j;
  if (!v)
    return;
  for (j = 0; j < n; j++)
    fq_clear(v + j, field);
  free(v);
}
int b_theta4_init(b_theta4_ctx *T, unsigned g, const fq_struct *null,
                  const fq_ctx_t field) {
  size_t j;
  memset(T, 0, sizeof(*T));
  if (!g || g > 8)
    return 0;
  T->g = g;
  T->count = (size_t)1 << (2 * g);
  T->characters = (size_t)1 << g;
  T->field = field;
  T->null = b_theta4_vector_new(T->count, field);
  T->charts = (b_theta4_chart *)calloc(T->count, sizeof(*T->charts));
  if (!T->null || !T->charts) {
    b_theta4_clear(T);
    return 0;
  }
  for (j = 0; j < T->count; j++) {
    fq_set(T->null + j, null + j, field);
    fq_init(T->charts[j].inverse_denominator, field);
  }
  return 1;
}
void b_theta4_clear(b_theta4_ctx *T) {
  size_t j;
  if (T->charts) {
    for (j = 0; j < T->count; j++)
      fq_clear(T->charts[j].inverse_denominator, T->field);
    free(T->charts);
  }
  b_theta4_vector_free(T->null, T->count, T->field);
  memset(T, 0, sizeof(*T));
}
void b_theta4_copy(fq_struct *out, const fq_struct *in, const b_theta4_ctx *T) {
  size_t j;
  for (j = 0; j < T->count; j++)
    fq_set(out + j, in + j, T->field);
}
int b_theta4_nonzero(const fq_struct *v, const b_theta4_ctx *T) {
  size_t j;
  for (j = 0; j < T->count; j++)
    if (!fq_is_zero(v + j, T->field))
      return 1;
  return 0;
}
int b_theta4_equal_projective(const fq_struct *a, const fq_struct *b,
                              const b_theta4_ctx *T) {
  size_t j, k;
  fq_t x, y;
  int ok = 1;
  for (k = 0; k < T->count; k++)
    if (!fq_is_zero(a + k, T->field))
      break;
  if (k == T->count || fq_is_zero(b + k, T->field))
    return 0;
  fq_init(x, T->field);
  fq_init(y, T->field);
  for (j = 0; j < T->count; j++) {
    fq_mul(x, a + j, b + k, T->field);
    fq_mul(y, b + j, a + k, T->field);
    if (!fq_equal(x, y, T->field)) {
      ok = 0;
      break;
    }
  }
  fq_clear(x, T->field);
  fq_clear(y, T->field);
  return ok;
}
void b_theta4_inverse(fq_struct *out, const fq_struct *in,
                      const b_theta4_ctx *T) {
  size_t j;
  fq_struct *v = b_theta4_vector_new(T->count, T->field);
  for (j = 0; j < T->count; j++)
    fq_set(v + j, in + negidx((uint32_t)j, T->g), T->field);
  b_theta4_copy(out, v, T);
  b_theta4_vector_free(v, T->count, T->field);
}
void b_theta4_action(fq_struct *out, const fq_struct *in, uint32_t alpha,
                     uint32_t beta, const fq_t zeta, const b_theta4_ctx *T) {
  size_t j;
  unsigned s, e;
  fq_t w;
  fq_struct *v = b_theta4_vector_new(T->count, T->field);
  fq_init(w, T->field);
  for (j = 0; j < T->count; j++) {
    e = 0;
    for (s = 0; s < T->g; s++)
      e += ((beta >> (2 * s)) & 3u) * (((uint32_t)j >> (2 * s)) & 3u);
    fq_pow_ui(w, zeta, e & 3u, T->field);
    fq_mul(v + j, w, in + addidx((uint32_t)j, alpha, T->g), T->field);
  }
  b_theta4_copy(out, v, T);
  b_theta4_vector_free(v, T->count, T->field);
  fq_clear(w, T->field);
}
static void L(fq_t out, unsigned chi, uint32_t i, uint32_t j,
              const fq_struct *P, const fq_struct *Q, const b_theta4_ctx *T) {
  size_t t;
  fq_t v;
  fq_init(v, T->field);
  fq_zero(out, T->field);
  for (t = 0; t < T->characters; t++) {
    uint32_t e = 2 * emb((uint32_t)t, T->g);
    fq_mul(v, P + addidx(i, e, T->g), Q + addidx(j, e, T->g), T->field);
    if (parity(chi & (unsigned)t))
      fq_sub(out, out, v, T->field);
    else
      fq_add(out, out, v, T->field);
  }
  fq_clear(v, T->field);
}
static void walsh(fq_struct *v,size_t n,const fq_ctx_t f);
static b_theta4_chart *chart_pair(unsigned chi, unsigned delta,
                                  b_theta4_chart *charts, const fq_struct *A,
                                  const fq_struct *B, b_theta4_ctx *T) {
  b_theta4_chart *c = charts + chi * T->characters + delta;
  uint32_t k, l;
  size_t s;
  if (c->ready)
    return c->ready > 0 ? c : NULL;
  fq_struct *den=b_theta4_vector_new(T->characters,T->field);
  if(!den)return NULL; /* Initial l has the required parity; all remaining
                             pairs follow. */
  for (k = 0; k < T->count; k++)
    for (s = 0; s < T->characters; s++) {
      l = emb(delta ^ lowbits(k, T->g), T->g) + 2 * emb((uint32_t)s, T->g);
      for(unsigned t=0;t<T->characters;t++){
        uint32_t e=2*emb(t,T->g);
        fq_mul(den+t,A+addidx(k,e,T->g),B+addidx(l,e,T->g),T->field);
      }
      walsh(den,T->characters,T->field);
      for(unsigned t=0;t<T->characters;t++){
        b_theta4_chart *d=charts+t*T->characters+delta;
        if(!d->ready&&!fq_is_zero(den+t,T->field)){
          d->k=k;d->l=l;fq_inv(d->inverse_denominator,den+t,T->field);d->ready=1;
        }
      }
      if(c->ready){b_theta4_vector_free(den,T->characters,T->field);return c;}
    }
  c->ready = -1;
  b_theta4_vector_free(den,T->characters,T->field);
  return NULL;
}
static b_theta4_chart *chart(unsigned chi, unsigned delta, b_theta4_ctx *T) {
  return chart_pair(chi, delta, T->charts, T->null, T->null, T);
}
static int riemann(fq_t out, unsigned chi, uint32_t i, uint32_t j,
                   const fq_struct *P, const fq_struct *Q, b_theta4_ctx *T) {
  b_theta4_chart *c = chart(chi, lowbits(i, T->g) ^ lowbits(j, T->g), T);
  uint32_t h, mh;
  fq_t a, b;
  if (!c)
    return 0;
  h = halfidx(i, j, c->k, c->l, T->g);
  mh = negidx(h, T->g);
  fq_init(a, T->field);
  fq_init(b, T->field);
  L(a, chi, addidx(i, h, T->g), addidx(j, mh, T->g), P, P, T);
  L(b, chi, addidx(c->k, mh, T->g), addidx(c->l, mh, T->g), Q, Q, T);
  fq_mul(out, a, b, T->field);
  fq_mul(out, out, c->inverse_denominator, T->field);
  fq_clear(a, T->field);
  fq_clear(b, T->field);
  return 1;
}
static void walsh(fq_struct *v, size_t n, const fq_ctx_t f) {
  size_t len, i, j;
  fq_t a, b;
  fq_init(a, f);
  fq_init(b, f);
  for (len = 1; len < n; len *= 2)
    for (i = 0; i < n; i += 2 * len)
      for (j = 0; j < len; j++) {
        fq_set(a, v + i + j, f);
        fq_set(b, v + i + j + len, f);
        fq_add(v + i + j, a, b, f);
        fq_sub(v + i + j + len, a, b, f);
      }
  fq_clear(a, f);
  fq_clear(b, f);
}
static int riemann_all(fq_struct *out,uint32_t i,uint32_t j,
    const fq_struct *P1,const fq_struct *P2,const fq_struct *Q1,const fq_struct *Q2,
    const fq_struct *denA,const fq_struct *denB,b_theta4_chart *charts,
    fq_struct *a,fq_struct *b,b_theta4_ctx *T){
  unsigned char used[256]={0};unsigned delta=lowbits(i,T->g)^lowbits(j,T->g);
  for(unsigned chi=0;chi<T->characters;chi++)if(!used[chi]){
    b_theta4_chart *c=chart_pair(chi,delta,charts,denA,denB,T);if(!c)return 0;
    uint32_t h=halfidx(i,j,c->k,c->l,T->g),mh=negidx(h,T->g);
    uint32_t ai=addidx(i,h,T->g),aj=addidx(j,mh,T->g);
    uint32_t bi=addidx(c->k,mh,T->g),bj=addidx(c->l,mh,T->g);
    for(unsigned s=0;s<T->characters;s++){
      uint32_t e=2*emb(s,T->g);
      fq_mul(a+s,P1+addidx(ai,e,T->g),P2+addidx(aj,e,T->g),T->field);
      fq_mul(b+s,Q1+addidx(bi,e,T->g),Q2+addidx(bj,e,T->g),T->field);
    }
    walsh(a,T->characters,T->field);walsh(b,T->characters,T->field);
    for(unsigned k=chi;k<T->characters;k++)if(!used[k]){
      b_theta4_chart *d=chart_pair(k,delta,charts,denA,denB,T);if(!d)return 0;
      if(d->k!=c->k||d->l!=c->l)continue;
      fq_mul(out+k,a+k,b+k,T->field);fq_mul(out+k,out+k,d->inverse_denominator,T->field);used[k]=1;
    }
  }
  return 1;
}
/* Compute the inverse-character sum without evaluating each L_chi
 * separately. Characters sharing a denominator (k,l) also share h. For
 * that group two Walsh transforms give all the required L values. The
 * finite chart fallback remains unchanged if groups split at exceptions. */
static int riemann_sum(fq_t out,uint32_t i,uint32_t j,
    const fq_struct *P,const fq_struct *Q,fq_struct *a,fq_struct *b,
    b_theta4_ctx *T){
  unsigned char used[256]={0};unsigned delta=lowbits(i,T->g)^lowbits(j,T->g);
  fq_t term;fq_init(term,T->field);fq_zero(out,T->field);int ok=0;
  for(unsigned chi=0;chi<T->characters;chi++)if(!used[chi]){
    b_theta4_chart *c=chart(chi,delta,T);if(!c)goto done;
    uint32_t h=halfidx(i,j,c->k,c->l,T->g),mh=negidx(h,T->g);
    uint32_t ai=addidx(i,h,T->g),aj=addidx(j,mh,T->g);
    uint32_t bi=addidx(c->k,mh,T->g),bj=addidx(c->l,mh,T->g);
    for(unsigned s=0;s<T->characters;s++){
      uint32_t e=2*emb(s,T->g);
      fq_mul(a+s,P+addidx(ai,e,T->g),P+addidx(aj,e,T->g),T->field);
      fq_mul(b+s,Q+addidx(bi,e,T->g),Q+addidx(bj,e,T->g),T->field);
    }
    walsh(a,T->characters,T->field);walsh(b,T->characters,T->field);
    for(unsigned k=chi;k<T->characters;k++)if(!used[k]){
      b_theta4_chart *d=chart(k,delta,T);if(!d)goto done;
      if(d->k!=c->k||d->l!=c->l)continue;
      fq_mul(term,a+k,b+k,T->field);fq_mul(term,term,d->inverse_denominator,T->field);
      fq_add(out,out,term,T->field);used[k]=1;
    }
  }
  ok=1;
done:fq_clear(term,T->field);return ok;
}
int b_theta4_add(fq_struct *out, const fq_struct *P, const fq_struct *Q,
                 b_theta4_ctx *T) {
  fq_struct *column=b_theta4_vector_new(T->count,T->field);
  fq_struct *a=b_theta4_vector_new(T->characters,T->field),*b=b_theta4_vector_new(T->characters,T->field);
  int ok=0;if(!column||!a||!b)goto done;
  for(size_t j=0;j<T->count;j++){
    for(size_t i=0;i<T->count;i++){
      if(!riemann_sum(column+i,(uint32_t)i,(uint32_t)j,P,Q,a,b,T))goto done;
    }
    /* This entire column is 2^g times (P-Q)_j times one lift of P+Q.
     * Keeping a fixed j is essential: diagonal products have different
     * factors and cannot be assembled into a single projective point. */
    if(b_theta4_nonzero(column,T)){
      b_theta4_copy(out,column,T);ok=1;break;
    }
  }
done:b_theta4_vector_free(column,T->count,T->field);b_theta4_vector_free(a,T->characters,T->field);
  b_theta4_vector_free(b,T->characters,T->field);return ok;
}
#ifndef B_THETA4_DIFF_IMPL
#define B_THETA4_DIFF_IMPL b_theta4_diff_add
#endif
int B_THETA4_DIFF_IMPL(fq_struct *out, const fq_struct *P, const fq_struct *Q,
                      const fq_struct *PmQ, b_theta4_ctx *T) {
  size_t q, s, i, pivot;
  int ok = 0;
  fq_t den, v, sum;
  fq_struct *tmp = b_theta4_vector_new(T->count, T->field),
            *block = b_theta4_vector_new(T->characters, T->field),
            *scratch = b_theta4_vector_new(2*T->characters,T->field);
  if (!tmp || !block || !scratch)
    goto allocfail;
  fq_init(den, T->field);
  fq_init(v, T->field);
  fq_init(sum, T->field);
  for (pivot = 0; pivot < T->count; pivot++)
    if (!fq_is_zero(PmQ + pivot, T->field))
      break;
  if (pivot == T->count)
    goto done;
  /* A single Walsh block supplies every simultaneous two-torsion shift. */
  for (q = 0; q < T->characters; q++) {
    uint32_t base = emb((uint32_t)q, T->g);
    if(!riemann_all(block,base,base,P,P,Q,Q,T->null,T->null,T->charts,
      scratch,scratch+T->characters,T))goto done;
    walsh(block, T->characters, T->field);
    for (s = 0; s < T->characters; s++) {
      i = base + 2 * emb((uint32_t)s, T->g);
      if (fq_is_zero(PmQ + i, T->field))
        continue;
      fq_mul_ui(den, PmQ + i, (ulong)T->characters, T->field);
      fq_inv(den, den, T->field);
      fq_mul(tmp + i, block + s, den, T->field);
    }
  }
  /* Coordinates where the difference vanishes use its first nonzero chart. */
  for (i = 0; i < T->count; i++)
    if (fq_is_zero(PmQ + i, T->field)) {
      if(!riemann_sum(sum,(uint32_t)i,(uint32_t)pivot,P,Q,
          scratch,scratch+T->characters,T))goto done;
      fq_mul_ui(den, PmQ + pivot, (ulong)T->characters, T->field);
      fq_inv(den, den, T->field);
      fq_mul(tmp + i, sum, den, T->field);
    }
  if (!b_theta4_nonzero(tmp, T))
    goto done;
  b_theta4_copy(out, tmp, T);
  ok = 1;
done:
  fq_clear(den, T->field);
  fq_clear(v, T->field);
  fq_clear(sum, T->field);
allocfail:
  b_theta4_vector_free(tmp, T->count, T->field);
  b_theta4_vector_free(block, T->characters, T->field);
  b_theta4_vector_free(scratch,2*T->characters,T->field);
  return ok;
}
/* Riemann position (P+Q+R,-P,-Q,-R ; 0,Q+R,P+R,P+Q). */
static int three_value(fq_t out, unsigned chi, uint32_t i, uint32_t j,
                       const fq_struct *nQ, const fq_struct *nR,
                       const fq_struct *QR, const fq_struct *PR,
                       const fq_struct *PQ, b_theta4_chart *charts,
                       b_theta4_ctx *T) {
  b_theta4_chart *c =
      chart_pair(chi, lowbits(i, T->g) ^ lowbits(j, T->g), charts, nQ, nR, T);
  uint32_t h, mh;
  fq_t a, b;
  if (!c)
    return 0;
  h = halfidx(i, j, c->k, c->l, T->g);
  mh = negidx(h, T->g);
  fq_init(a, T->field);
  fq_init(b, T->field);
  L(a, chi, addidx(i, h, T->g), addidx(j, mh, T->g), T->null, QR, T);
  L(b, chi, addidx(c->k, mh, T->g), addidx(c->l, mh, T->g), PR, PQ, T);
  fq_mul(out, a, b, T->field);
  fq_mul(out, out, c->inverse_denominator, T->field);
  fq_clear(a, T->field);
  fq_clear(b, T->field);
  return 1;
}
#ifndef B_THETA4_THREE_IMPL
#define B_THETA4_THREE_IMPL b_theta4_three_way_add
#endif
int B_THETA4_THREE_IMPL(fq_struct *out, const fq_struct *P,
                           const fq_struct *Q, const fq_struct *R,
                           const fq_struct *PQ, const fq_struct *QR,
                           const fq_struct *PR, b_theta4_ctx *T) {
  size_t i, q, chi, s, pivot;
  int ok = 0;
  fq_t den, v, sum;
  fq_struct *nP = b_theta4_vector_new(T->count, T->field),
            *nQ = b_theta4_vector_new(T->count, T->field),
            *nR = b_theta4_vector_new(T->count, T->field),
            *tmp = b_theta4_vector_new(T->count, T->field),
            *block = b_theta4_vector_new(T->characters, T->field),
            *scratch = b_theta4_vector_new(2*T->characters,T->field);
  b_theta4_chart *charts = (b_theta4_chart *)calloc(T->count, sizeof(*charts));
  if (!nP || !nQ || !nR || !tmp || !block || !scratch || !charts)
    goto allocfail;
  for (i = 0; i < T->count; i++)
    fq_init(charts[i].inverse_denominator, T->field);
  fq_init(den, T->field);
  fq_init(v, T->field);
  fq_init(sum, T->field);
  b_theta4_inverse(nP, P, T);
  b_theta4_inverse(nQ, Q, T);
  b_theta4_inverse(nR, R, T);
  for (pivot = 0; pivot < T->count; pivot++)
    if (!fq_is_zero(nP + pivot, T->field))
      break;
  if (pivot == T->count)
    goto done;
  for (q = 0; q < T->characters; q++) {
    uint32_t base = emb((uint32_t)q, T->g);
    if(!riemann_all(block,base,base,T->null,QR,PR,PQ,nQ,nR,charts,
      scratch,scratch+T->characters,T))goto done;
    walsh(block, T->characters, T->field);
    for (s = 0; s < T->characters; s++) {
      i = base + 2 * emb((uint32_t)s, T->g);
      if (fq_is_zero(nP + i, T->field))
        continue;
      fq_mul_ui(den, nP + i, (ulong)T->characters, T->field);
      fq_inv(den, den, T->field);
      fq_mul(tmp + i, block + s, den, T->field);
    }
  }
  for (i = 0; i < T->count; i++)
    if (fq_is_zero(nP + i, T->field)) {
      if(!riemann_all(block,(uint32_t)i,(uint32_t)pivot,T->null,QR,PR,PQ,
          nQ,nR,charts,scratch,scratch+T->characters,T))goto done;
      fq_zero(sum,T->field);
      for(chi=0;chi<T->characters;chi++)fq_add(sum,sum,block+chi,T->field);
      fq_mul_ui(den, nP + pivot, (ulong)T->characters, T->field);
      fq_inv(den, den, T->field);
      fq_mul(tmp + i, sum, den, T->field);
    }
  if (!b_theta4_nonzero(tmp, T))
    goto done;
  b_theta4_copy(out, tmp, T);
  ok = 1;
done:
  fq_clear(den, T->field);
  fq_clear(v, T->field);
  fq_clear(sum, T->field);
  for (i = 0; i < T->count; i++)
    fq_clear(charts[i].inverse_denominator, T->field);
allocfail:
  free(charts);
  b_theta4_vector_free(nP, T->count, T->field);
  b_theta4_vector_free(nQ, T->count, T->field);
  b_theta4_vector_free(nR, T->count, T->field);
  b_theta4_vector_free(tmp, T->count, T->field);
  b_theta4_vector_free(block, T->characters, T->field);
  b_theta4_vector_free(scratch,2*T->characters,T->field);
  return ok;
}
int b_theta4_scalar(fq_struct *out, const fq_struct *P, ulong n,
                    b_theta4_ctx *T) {
  fq_struct *a = b_theta4_vector_new(T->count, T->field),
            *b = b_theta4_vector_new(T->count, T->field),
            *c = b_theta4_vector_new(T->count, T->field),
            *d = b_theta4_vector_new(T->count, T->field);
  ulong bit = 1;
  int ok = 0;
  if (!a || !b || !c || !d)
    goto done;
  b_theta4_copy(a, T->null, T);
  b_theta4_copy(b, P, T);
  while (bit <= n / 2)
    bit *= 2;
  for (; bit; bit /= 2) {
    if (!b_theta4_diff_add(c, b, a, P, T))
      goto done;
    if (n & bit) {
      if (!b_theta4_diff_add(d, b, b, T->null, T))
        goto done;
      b_theta4_copy(a, c, T);
      b_theta4_copy(b, d, T);
    } else {
      if (!b_theta4_diff_add(d, a, a, T->null, T))
        goto done;
      b_theta4_copy(a, d, T);
      b_theta4_copy(b, c, T);
    }
  }
  b_theta4_copy(out, a, T);
  ok = 1;
done:
  b_theta4_vector_free(a, T->count, T->field);
  b_theta4_vector_free(b, T->count, T->field);
  b_theta4_vector_free(c, T->count, T->field);
  b_theta4_vector_free(d, T->count, T->field);
  return ok;
}
int b_theta4_scalar_add(fq_struct *nPQ, fq_struct *nP, ulong n,
                        const fq_struct *PQ, const fq_struct *P,
                        const fq_struct *Q, b_theta4_ctx *T) {
  fq_struct *v[8] = {0};
  ulong bit = 1;
  int ok = 0;
  unsigned i;
  for (i = 0; i < 8; i++)
    if (!(v[i] = b_theta4_vector_new(T->count, T->field)))
      goto done;
  /* v0=kP,v1=(k+1)P,v2=kP+Q,v3=(k+1)P+Q. */
  b_theta4_copy(v[0], T->null, T);
  b_theta4_copy(v[1], P, T);
  b_theta4_copy(v[2], Q, T);
  b_theta4_copy(v[3], PQ, T);
  while (bit <= n / 2)
    bit *= 2;
  for (; bit; bit /= 2) {
    if (!b_theta4_diff_add(v[4], v[1], v[0], P, T) ||
        !b_theta4_diff_add(v[6], v[3], v[0], PQ, T))
      goto done;
    if (n & bit) {
      if (!b_theta4_diff_add(v[5], v[1], v[1], T->null, T) ||
          !b_theta4_diff_add(v[7], v[3], v[1], Q, T))
        goto done;
      b_theta4_copy(v[0], v[4], T);
      b_theta4_copy(v[1], v[5], T);
      b_theta4_copy(v[2], v[6], T);
      b_theta4_copy(v[3], v[7], T);
    } else {
      if (!b_theta4_diff_add(v[5], v[0], v[0], T->null, T) ||
          !b_theta4_diff_add(v[7], v[2], v[0], Q, T))
        goto done;
      b_theta4_copy(v[0], v[5], T);
      b_theta4_copy(v[1], v[4], T);
      b_theta4_copy(v[2], v[7], T);
      b_theta4_copy(v[3], v[6], T);
    }
  }
  b_theta4_copy(nPQ, v[2], T);
  b_theta4_copy(nP, v[0], T);
  ok = 1;
done:
  for (i = 0; i < 8; i++)
    b_theta4_vector_free(v[i], T->count, T->field);
  return ok;
}
int b_theta4_scale_factor(fq_t c, const fq_struct *source,
                          const fq_struct *target, const b_theta4_ctx *T) {
  size_t i, k;
  fq_t scale, w;
  int ok = 0;
  for (k = 0; k < T->count; k++)
    if (!fq_is_zero(source + k, T->field))
      break;
  if (k == T->count || fq_is_zero(target + k, T->field))
    return 0;
  fq_init(scale, T->field);
  fq_init(w, T->field);
  fq_inv(scale, source + k, T->field);
  fq_mul(scale, scale, target + k, T->field);
  for (i = 0; i < T->count; i++) {
    fq_mul(w, source + i, scale, T->field);
    if (!fq_equal(w, target + i, T->field))
      goto done;
  }
  fq_set(c, scale, T->field);
  ok = 1;
done:
  fq_clear(scale, T->field);
  fq_clear(w, T->field);
  return ok;
}
int b_theta4_good_basis2(fq_struct *out, const fq_struct *P, uint32_t alpha,
                         const fq_t zeta, b_theta4_ctx *T) {
  fq_struct *v[3] = {0};
  fq_t scale, s, lambda;
  unsigned i;
  int status = 0;
  for (i = 0; i < 3; i++)
    if (!(v[i] = b_theta4_vector_new(T->count, T->field)))
      goto allocfail;
  fq_init(scale, T->field);
  fq_init(s, T->field);
  fq_init(lambda, T->field);
  b_theta4_inverse(v[0], P, T);
  b_theta4_action(v[1], v[0], alpha, 0, zeta, T);
  if (!b_theta4_scale_factor(scale, P, v[1], T))
    goto done;
  if (!fq_is_one(scale, T->field)) {
    status = -1;
    goto done;
  }
  if (!b_theta4_diff_add(v[0], P, P, T->null, T))
    goto done;
  b_theta4_action(v[1], T->null, alpha, 0, zeta, T);
  if (!b_theta4_scale_factor(scale, v[0], v[1], T))
    goto done;
  if (!fq_sqrt(s, scale, T->field)) {
    status = -2;
    goto done;
  }
  if (!fq_sqrt(lambda, s, T->field)) {
    fq_neg(s, s, T->field);
    if (!fq_sqrt(lambda, s, T->field)) {
      status = -2;
      goto done;
    }
  }
  for (i = 0; i < T->count; i++)
    fq_mul(v[2] + i, P + i, lambda, T->field);
  if (!b_theta4_diff_add(v[0], v[2], v[2], T->null, T))
    goto done;
  for (i = 0; i < T->count; i++)
    if (!fq_equal(v[0] + i, v[1] + i, T->field))
      goto done;
  b_theta4_copy(out, v[2], T);
  status = 1;
done:
  fq_clear(scale, T->field);
  fq_clear(s, T->field);
  fq_clear(lambda, T->field);
allocfail:
  for (i = 0; i < 3; i++)
    b_theta4_vector_free(v[i], T->count, T->field);
  return status;
}
int b_theta4_good_coset2(fq_struct *out, const fq_struct *PQ,
                         const fq_struct *P, const fq_struct *Q, uint32_t alpha,
                         const fq_t zeta, b_theta4_ctx *T) {
  fq_struct *v[3] = {0};
  fq_t scale, lambda;
  unsigned i;
  int status = 0;
  for (i = 0; i < 3; i++)
    if (!(v[i] = b_theta4_vector_new(T->count, T->field)))
      goto allocfail;
  fq_init(scale, T->field);
  fq_init(lambda, T->field);
  if (!b_theta4_diff_add(v[0], PQ, P, Q, T))
    goto done;
  b_theta4_action(v[1], Q, alpha, 0, zeta, T);
  if (!b_theta4_scale_factor(scale, v[0], v[1], T))
    goto done;
  if (!fq_sqrt(lambda, scale, T->field)) {
    status = -2;
    goto done;
  }
  for (i = 0; i < T->count; i++)
    fq_mul(v[2] + i, PQ + i, lambda, T->field);
  if (!b_theta4_diff_add(v[0], v[2], P, Q, T))
    goto done;
  for (i = 0; i < T->count; i++)
    if (!fq_equal(v[0] + i, v[1] + i, T->field))
      goto done;
  b_theta4_copy(out, v[2], T);
  status = 1;
done:
  fq_clear(scale, T->field);
  fq_clear(lambda, T->field);
allocfail:
  for (i = 0; i < 3; i++)
    b_theta4_vector_free(v[i], T->count, T->field);
  return status;
}
int b_theta4_good_coset2_squared_scale(fq_t s,const fq_struct *PQ,
    const fq_struct *P,const fq_t sP,const fq_struct *Q,uint32_t alpha,
    const fq_t zeta,b_theta4_ctx *T) {
  fq_struct *v[2]={0};fq_t scale;int ok=0;
  if(fq_is_zero(sP,T->field))return 0;
  v[0]=b_theta4_vector_new(T->count,T->field);
  v[1]=b_theta4_vector_new(T->count,T->field);
  if(!v[0]||!v[1])goto done;
  fq_init(scale,T->field);
  if(b_theta4_diff_add(v[0],PQ,P,Q,T)){
    b_theta4_action(v[1],Q,alpha,0,zeta,T);
    if(b_theta4_scale_factor(scale,v[0],v[1],T)){
      fq_div(scale,scale,sP,T->field);fq_set(s,scale,T->field);ok=1;
    }
  }
  fq_clear(scale,T->field);
done:
  b_theta4_vector_free(v[0],T->count,T->field);
  b_theta4_vector_free(v[1],T->count,T->field);return ok;
}
void b_theta4_two_standard_block(fq_struct *out,unsigned parity,
    const fq_struct *shifted,const b_theta4_ctx *T) {
  size_t j,k;fq_t sum,square;fq_init(sum,T->field);fq_init(square,T->field);
  for(j=0;j<T->characters;j++){
    uint32_t jj=emb((uint32_t)j,T->g);
    fq_zero(sum,T->field);
    for(k=0;k<T->characters;k++){
      fq_sqr(square,shifted+jj+2*emb((uint32_t)k,T->g),T->field);
      fq_add(sum,sum,square,T->field);
    }
    fq_set(out+emb(parity,T->g)+2*jj,sum,T->field);
  }
  fq_clear(square,T->field);fq_clear(sum,T->field);
}
void b_theta4_two_standard_scaled_block(fq_struct *out,unsigned parity,
    const fq_struct *shifted,const fq_t squared_scale,const b_theta4_ctx *T){
  b_theta4_two_standard_block(out,parity,shifted,T);
  for(size_t j=0;j<T->characters;j++){
    size_t index=emb(parity,T->g)+2*emb((uint32_t)j,T->g);
    fq_mul(out+index,out+index,squared_scale,T->field);
  }
}
void b_theta4_two_trace_block(fq_struct *out,unsigned parity,unsigned seed,
    const fq_struct *shifted,const fq_t squared_scale,const b_theta4_ctx *T){
  if(!seed){b_theta4_two_standard_scaled_block(out,parity,shifted,squared_scale,T);return;}
  uint32_t plus=emb(seed,T->g),minus=negidx(plus,T->g);
  fq_t sum,term;fq_init(sum,T->field);fq_init(term,T->field);
  for(unsigned j=0;j<T->characters;j++){
    uint32_t q=emb(j,T->g);fq_zero(sum,T->field);
    for(unsigned k=0;k<T->characters;k++){
      uint32_t index=q+2*emb(k,T->g);
      fq_mul(term,shifted+addidx(index,plus,T->g),shifted+addidx(index,minus,T->g),T->field);
      fq_add(sum,sum,term,T->field);
    }
    fq_mul(out+emb(parity,T->g)+2*q,sum,squared_scale,T->field);
  }
  fq_clear(sum,T->field);fq_clear(term,T->field);
}
int b_theta4_to_kummer(fq_t U, fq_t V, const fq_struct *P,
                       const b_theta4_ctx *T) {
  fq_t a, b, c, w;
  int ok = 0;
  if (T->g != 1)
    return 0;
  fq_init(a, T->field);
  fq_init(b, T->field);
  fq_init(c, T->field);
  fq_init(w, T->field);
  fq_mul(a, P, T->null, T->field);
  fq_mul(w, P + 2, T->null + 2, T->field);
  fq_add(a, a, w, T->field);
  fq_mul(b, P + 1, T->null + 1, T->field);
  fq_mul(w, P + 3, T->null + 3, T->field);
  fq_add(b, b, w, T->field);
  fq_mul(c, P, T->null + 2, T->field);
  fq_mul(w, P + 2, T->null, T->field);
  fq_add(c, c, w, T->field);
  if (!fq_is_zero(a, T->field)) {
    fq_set(U, a, T->field);
    fq_set(V, b, T->field);
    ok = 1;
  } else if (!fq_is_zero(c, T->field)) {
    fq_set(U, b, T->field);
    fq_set(V, c, T->field);
    ok = 1;
  }
  fq_clear(a, T->field);
  fq_clear(b, T->field);
  fq_clear(c, T->field);
  fq_clear(w, T->field);
  return ok;
}
void b_theta4_fourier_axis(fq_struct *out, const fq_struct *in, unsigned axis,
                           const fq_t zeta, const b_theta4_ctx *T) {
  size_t i;
  unsigned j, digit;
  uint32_t mask;
  fq_t w, v;
  fq_struct *tmp = b_theta4_vector_new(T->count, T->field);
  if (axis >= T->g || !tmp)
    return;
  mask = 3u << (2 * axis);
  fq_init(w, T->field);
  fq_init(v, T->field);
  for (i = 0; i < T->count; i++) {
    digit = ((uint32_t)i >> (2 * axis)) & 3u;
    fq_zero(tmp + i, T->field);
    for (j = 0; j < 4; j++) {
      fq_pow_ui(w, zeta, (digit * j) & 3u, T->field);
      fq_mul(v, w, in + (((uint32_t)i & ~mask) | (j << (2 * axis))), T->field);
      fq_add(tmp + i, tmp + i, v, T->field);
    }
  }
  b_theta4_copy(out, tmp, T);
  b_theta4_vector_free(tmp, T->count, T->field);
  fq_clear(w, T->field);
  fq_clear(v, T->field);
}
