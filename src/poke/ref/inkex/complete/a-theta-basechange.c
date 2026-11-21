#include "a-theta-basechange.h"
#include <stdlib.h>
#include <string.h>

void a_theta_two_basechange_clear(a_theta_two_basechange *S)
{
    if(!S||!S->ready)return;
    b_theta4_two_clear(&S->quotient);b_theta4_clear(&S->source);
    memset(S,0,sizeof(*S));
}

int a_theta_two_basechange_init(a_theta_two_basechange *out,
    const b_theta4_two_ctx *Q,const root_field_embedding *L)
{
    if(!out||out->ready||!Q||!Q->source||!Q->kernel_lift||!Q->target_null||
       !L||!L->ready||Q->source->field!=L->sub->field)return 0;
    size_t n=Q->source->count;const fq_ctx_struct *field=L->sup->field;
    fq_struct *null=b_theta4_vector_new(n,field);if(!null)return 0;
    for(size_t i=0;i<n;i++)root_field_up(null+i,Q->source->null+i,L);
    a_theta_two_basechange tmp={0};int ok=0;
    if(!b_theta4_init(&tmp.source,Q->source->g,null,field))goto done;
    tmp.ready=1;tmp.quotient.source=&tmp.source;
    tmp.quotient.kernel_lift=calloc(tmp.source.characters,sizeof(*tmp.quotient.kernel_lift));
    tmp.quotient.target_null=b_theta4_vector_new(n,field);
    if(!tmp.quotient.kernel_lift||!tmp.quotient.target_null)goto done;
    for(size_t s=0;s<tmp.source.characters;s++)if(Q->kernel_lift[s]){
        tmp.quotient.kernel_lift[s]=b_theta4_vector_new(n,field);
        if(!tmp.quotient.kernel_lift[s])goto done;
        for(size_t i=0;i<n;i++)root_field_up(tmp.quotient.kernel_lift[s]+i,Q->kernel_lift[s]+i,L);
    }
    for(size_t i=0;i<n;i++)root_field_up(tmp.quotient.target_null+i,Q->target_null+i,L);
    tmp.quotient.trace_seed=Q->trace_seed;
    *out=tmp;out->quotient.source=&out->source;memset(&tmp,0,sizeof(tmp));ok=1;
done:a_theta_two_basechange_clear(&tmp);b_theta4_vector_free(null,n,field);return ok;
}
