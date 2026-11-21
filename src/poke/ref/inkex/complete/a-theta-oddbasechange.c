#include "a-theta-oddbasechange.h"
#include <string.h>
void a_theta_odd_basechange_clear(a_theta_odd_basechange *S){if(!S||!S->ready)return;b_theta4_odd_clear(&S->quotient);b_theta4_clear(&S->source);memset(S,0,sizeof(*S));}
int a_theta_odd_basechange_init(a_theta_odd_basechange *out,const b_theta4_odd_ctx *Q,const root_field_embedding *L){
 if(!out||out->ready||!Q||!Q->source||!Q->basis||!Q->target_null||!L||!L->ready||Q->source->field!=L->sub->field)return 0;
 size_t n=Q->source->count,g=Q->source->g;const fq_ctx_struct *F=L->sup->field;
 fq_struct *null=b_theta4_vector_new(n,F);if(!null)return 0;for(size_t i=0;i<n;i++)root_field_up(null+i,Q->source->null+i,L);
 a_theta_odd_basechange S={0};int ok=0;if(!b_theta4_init(&S.source,(unsigned)g,null,F))goto done;S.ready=1;
 b_theta4_odd_ctx *R=&S.quotient;R->source=&S.source;R->ell=Q->ell;R->kernel_size=Q->kernel_size;R->selected_chart=Q->selected_chart;
 R->basis=b_theta4_vector_new(g*n,F);R->translate=b_theta4_vector_new(n,F);R->weights=b_theta4_vector_new(g,F);R->target_null=b_theta4_vector_new(n,F);
 if(!R->basis||!R->translate||!R->weights||!R->target_null)goto done;
 for(size_t i=0;i<g*n;i++)root_field_up(R->basis+i,Q->basis+i,L);
 for(size_t i=0;i<g;i++)root_field_up(R->weights+i,Q->weights+i,L);
 for(size_t i=0;i<n;i++){root_field_up(R->translate+i,Q->translate+i,L);root_field_up(R->target_null+i,Q->target_null+i,L);}
 *out=S;out->quotient.source=&out->source;memset(&S,0,sizeof(S));ok=1;
done:a_theta_odd_basechange_clear(&S);b_theta4_vector_free(null,n,F);return ok;
}
