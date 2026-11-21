#include "b-theta-balanced.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

void b_theta_balanced_clear(b_theta_balanced *B){
    if(!B)return;
    a_theta_multifield_clear(&B->state);free(B->height);memset(B,0,sizeof(*B));
}
int b_theta_balanced_init(b_theta_balanced *B,unsigned g,unsigned exponent,
    const root_extension *base,const root_extension *kernel_field,
    const fq_struct *null,const fq_struct *const *generators,
    const a_theta_payload_input *payloads,size_t group_count,
    const fq_t zeta4,const fq_t zeta8){
    if(!B||B->ready||!exponent||!g||g>8||group_count==SIZE_MAX||
       group_count>SIZE_MAX/sizeof(a_theta_payload_input)-1)return 0;
    a_theta_payload_input *all=calloc(group_count+1,sizeof(*all));if(!all)return 0;
    if(group_count){if(!payloads){free(all);return 0;}memcpy(all,payloads,group_count*sizeof(*all));}
    all[group_count].field=kernel_field;all[group_count].points=generators;all[group_count].count=g;
    int ok=a_theta_multifield_init(&B->state,g,exponent,base,kernel_field,null,generators,all,group_count+1,zeta4,zeta8);free(all);if(!ok)return 0;
    unsigned capacity=1,t=exponent;while(t>1){t/=2;++capacity;}
    B->height=calloc(capacity,sizeof(*B->height));if(!B->height){a_theta_multifield_clear(&B->state);return 0;}
    B->stack_group=group_count;B->capacity=capacity;B->depth=B->max_depth=1;B->height[0]=exponent;B->ready=1;return 1;
}
static int push_frame(b_theta_balanced *B){
    a_theta_multifield *S=&B->state;a_theta_payload_group *G=S->groups+B->stack_group;
    unsigned g=S->theta.g,h=B->height[B->depth-1],child=h/2,shift=h-child;size_t n=S->theta.count;
    if(h<2||B->depth>=B->capacity||G->count>SIZE_MAX/n-g)return 0;
    const root_extension *E=S->kernel_field;size_t old_count=G->count;
    fq_struct *next=b_theta4_vector_new((old_count+g)*n,E->field),*null=b_theta4_vector_new(n,E->field);
    b_theta4_ctx T={0};int ready=0,ok=0;
    if(!next||!null)goto done;
    for(size_t i=0;i<n;i++)if(E==S->base)fq_set(null+i,S->theta.null+i,E->field);else root_field_up(null+i,S->theta.null+i,&S->kernel_embedding);
    if(!b_theta4_init(&T,g,null,E->field))goto done;ready=1;
    for(size_t i=0;i<old_count*n;i++)fq_set(next+i,G->points+i,E->field);
    for(unsigned j=0;j<g;j++){
        fq_struct *P=next+(old_count+j)*n;const fq_struct *parent=G->points+(old_count-g+j)*n;
        b_theta4_copy(P,parent,&T);
        for(unsigned k=0;k<shift;k++)if(!b_theta4_diff_add(P,P,P,T.null,&T))goto done;
    }
    b_theta4_vector_free(G->points,old_count*n,E->field);G->points=next;next=NULL;G->count+=g;
    B->height[B->depth++]=child;if(B->depth>B->max_depth)B->max_depth=B->depth;
    B->frame_doublings+=(unsigned long long)g*shift;ok=1;
done:if(ready)b_theta4_clear(&T);b_theta4_vector_free(null,n,E->field);b_theta4_vector_free(next,(old_count+g)*n,E->field);return ok;
}
int b_theta_balanced_step(b_theta_balanced *B){
    if(!B||!B->ready||!B->state.remaining||!B->depth)return 0;
    a_theta_multifield *S=&B->state;unsigned g=S->theta.g;size_t n=S->theta.count;
    while(B->height[B->depth-1]>1)if(!push_frame(B))return 0;
    a_theta_payload_group *G=S->groups+B->stack_group;
    /* Allocate the compacted stack before the quotient, preserving the
     * external atomic-stage guarantee if allocation fails. */
    size_t after_count=G->count-g;
    fq_struct *after=after_count?b_theta4_vector_new(after_count*n,S->kernel_field->field):NULL;
    if(after_count&&!after)return 0;
    for(size_t i=0;i<(size_t)g*n;i++)fq_set(S->generators+i,G->points+after_count*n+i,S->kernel_field->field);
    unsigned remaining=S->remaining;S->remaining=1;
    int ok=a_theta_multifield_step(S);S->remaining=ok?remaining-1:remaining;
    if(!ok){b_theta4_vector_free(after,after_count*n,S->kernel_field->field);return 0;}
    G=S->groups+B->stack_group;
    for(size_t i=0;i<after_count*n;i++)fq_set(after+i,G->points+i,S->kernel_field->field);
    b_theta4_vector_free(G->points,G->count*n,S->kernel_field->field);G->points=after;G->count=after_count;
    B->frame_evaluations+=(unsigned long long)(B->depth+1)*g; /* stack plus current order-8 basis */
    for(unsigned i=0;i<B->depth;i++)B->height[i]--;
    B->depth--;return 1;
}
