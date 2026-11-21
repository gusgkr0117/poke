#include "a-theta-multifield.h"
#include "root-theta-label.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

void a_theta_multifield_clear(a_theta_multifield *S)
{
    if(!S||!S->ready)return;
    size_t n=S->theta.count;
    b_theta4_vector_free(S->generators,S->theta.g*n,S->base->field);
    for(size_t j=0;j<S->group_count;j++){
        a_theta_payload_group *G=S->groups+j;
        b_theta4_vector_free(G->points,G->count*n,S->base->field);
    }
    free(S->groups);
    fq_clear(S->zeta8,S->base->field);fq_clear(S->zeta4,S->base->field);
    b_theta4_clear(&S->theta);memset(S,0,sizeof(*S));
}
int a_theta_multifield_init(a_theta_multifield *S,unsigned g,unsigned exponent,
    const root_extension *base,
    const fq_struct *null,const fq_struct *const *generators,
    const a_theta_payload_input *payloads,size_t group_count,
    const fq_t zeta4,const fq_t zeta8)
{
    if(!S||S->ready||!base||!g||g>8||!exponent||
       !null||!generators||(group_count&&!payloads)||
       fq_ctx_degree(base->field)!=2)return 0;
    fq_t w;fq_init(w,base->field);fq_sqr(w,zeta4,base->field);fq_neg(w,w,base->field);
    int valid=fq_is_one(w,base->field);fq_sqr(w,zeta8,base->field);
    valid=valid&&fq_equal(w,zeta4,base->field);fq_clear(w,base->field);
    if(!valid||!b_theta4_init(&S->theta,g,null,base->field))return 0;
    S->base=base;S->remaining=exponent;
    fq_init(S->zeta4,base->field);fq_init(S->zeta8,base->field);S->ready=1;
    fq_set(S->zeta4,zeta4,base->field);fq_set(S->zeta8,zeta8,base->field);
    size_t n=S->theta.count;S->generators=b_theta4_vector_new(g*n,base->field);
    if(!S->generators)goto fail;
    for(unsigned j=0;j<g;j++){
        if(!generators[j])goto fail;
        for(size_t i=0;i<n;i++)fq_set(S->generators+j*n+i,generators[j]+i,base->field);
    }
    if(group_count){S->groups=calloc(group_count,sizeof(*S->groups));if(!S->groups)goto fail;}
    S->group_count=group_count;
    for(size_t j=0;j<group_count;j++){
        const a_theta_payload_input *P=payloads+j;a_theta_payload_group *G=S->groups+j;
        if((P->count&&!P->points)||P->count>SIZE_MAX/n)goto fail;
        G->count=P->count;
        if(G->count){G->points=b_theta4_vector_new(G->count*n,S->base->field);if(!G->points)goto fail;}
        for(size_t k=0;k<G->count;k++){
            if(!P->points[k])goto fail;
            for(size_t i=0;i<n;i++)fq_set(G->points+k*n+i,P->points[k]+i,S->base->field);
        }
    }
    return 1;
fail:
    a_theta_multifield_clear(S);
    return 0;
}
static int mark_in_field(fq_struct *out,const fq_struct *in,const root_sp_word *word,
    unsigned mask,unsigned g,const fq_t zeta8,const root_extension *E)
{
    return root_sp_apply(out,in,word,zeta8,0,E->field)&&
        root_theta4_compatibility_apply(out,out,mask,g,E->field);
}
int a_theta_multifield_step(a_theta_multifield *S)
{
    if(!S||!S->ready||!S->remaining)return 0;
    const root_extension *B=S->base;unsigned g=S->theta.g;
    size_t n=S->theta.count;int ok=0,source_ready=0,target_ready=0;
    fq_struct *new_generators=b_theta4_vector_new(g*n,B->field),*raw=b_theta4_vector_new(n,B->field);
    fq_struct *null=b_theta4_vector_new(n,B->field),*basis=b_theta4_vector_new(g*n,B->field);
    fq_struct *raw_basis=b_theta4_vector_new(g*n,B->field),*pairs[64]={0};
    fq_struct **new_groups=S->group_count?calloc(S->group_count,sizeof(*new_groups)):NULL;
    b_theta4_ctx source={0},target={0};b_theta4_two_ctx quotient={0};
    root_theta_label_ctx labels={0};
    root_sp_matrix columns={{0}},frame,inverse;root_sp_word word;unsigned mask=0;
    if(!new_generators||!raw||!null||!basis||!raw_basis||
       (S->group_count&&!new_groups))goto done;
    if(!root_theta_label_init(&labels,&S->theta,S->zeta4))goto done;
    for(unsigned j=0;j<g;j++){
        fq_struct *generator=raw_basis+j*n;
        b_theta4_copy(generator,S->generators+j*n,&S->theta);
        for(unsigned k=1;k<S->remaining;k++)
            if(!b_theta4_diff_add(generator,generator,generator,S->theta.null,&S->theta))goto done;
        if(!b_theta4_diff_add(raw,generator,generator,S->theta.null,&S->theta))goto done;
        uint32_t alpha,beta;if(!root_theta_label(&alpha,&beta,raw,&labels))goto done;
        for(unsigned i=0;i<g;i++){columns[i][j]=(alpha>>(2*i))&3;columns[g+i][j]=(beta>>(2*i))&3;}
    }
    if(!root_sp_complete_lagrangian(frame,columns,g)||!root_sp_inverse(inverse,frame,g)||
       !root_sp_decompose(&word,inverse,g)||!mark_in_field(null,S->theta.null,&word,0,g,S->zeta8,B))goto done;
    const fq_struct *bases[8];for(unsigned j=0;j<g;j++){
        if(!mark_in_field(basis+j*n,raw_basis+j*n,&word,0,g,S->zeta8,B))goto done;bases[j]=basis+j*n;
        for(unsigned k=j+1;k<g;k++){
            pairs[j*g+k]=b_theta4_vector_new(n,B->field);
            if(!pairs[j*g+k]||!b_theta4_add(raw,raw_basis+j*n,raw_basis+k*n,&S->theta)||
               !mark_in_field(pairs[j*g+k],raw,&word,0,g,S->zeta8,B))goto done;
        }
    }
    if(!b_theta4_init(&source,g,null,B->field))goto done;source_ready=1;
    for(unsigned j=0;j<g;j++){int sign;if(!root_theta4_compatibility_sign(&sign,basis+j*n,1u<<(2*j),&source))goto done;if(sign<0)mask|=1u<<j;}
    if(!root_theta4_compatibility_apply(null,null,mask,g,B->field))goto done;
    for(unsigned j=0;j<g;j++){
        if(!root_theta4_compatibility_apply(basis+j*n,basis+j*n,mask,g,B->field))goto done;
        for(unsigned k=j+1;k<g;k++)if(!root_theta4_compatibility_apply(pairs[j*g+k],pairs[j*g+k],mask,g,B->field))goto done;
    }
    b_theta4_clear(&source);source_ready=0;if(!b_theta4_init(&source,g,null,B->field))goto done;source_ready=1;
    if(b_theta4_two_init(&quotient,&source,bases,(const fq_struct*const*)pairs,S->zeta4)!=1||
       !b_theta4_init(&target,g,quotient.target_null,B->field))goto done;
    target_ready=1;
    for(unsigned j=0;j<g;j++)
        if(!mark_in_field(raw,S->generators+j*n,&word,mask,g,S->zeta8,B)||
           !b_theta4_two_evaluate_point(new_generators+j*n,raw,S->zeta4,&quotient))goto done;
    for(size_t k=0;k<S->group_count;k++){
        a_theta_payload_group *G=S->groups+k;
        if(!G->count)continue;
        new_groups[k]=b_theta4_vector_new(G->count*n,B->field);
        if(!new_groups[k])goto done;
        for(size_t j=0;j<G->count;j++)
            if(!mark_in_field(raw,G->points+j*n,&word,mask,g,S->zeta8,B)||
               !b_theta4_two_evaluate_point(new_groups[k]+j*n,raw,S->zeta4,&quotient))goto done;
    }
    root_theta_label_clear(&labels);
    b_theta4_clear(&S->theta);S->theta=target;memset(&target,0,sizeof(target));target_ready=0;
    b_theta4_vector_free(S->generators,g*n,B->field);S->generators=new_generators;new_generators=NULL;
    for(size_t j=0;j<S->group_count;j++){
        a_theta_payload_group *G=S->groups+j;b_theta4_vector_free(G->points,G->count*n,S->base->field);
        G->points=new_groups[j];new_groups[j]=NULL;
    }
    S->remaining--;S->completed++;ok=1;
done:
    root_theta_label_clear(&labels);b_theta4_two_clear(&quotient);
    if(target_ready)b_theta4_clear(&target);if(source_ready)b_theta4_clear(&source);
    if(new_groups){for(size_t j=0;j<S->group_count;j++)b_theta4_vector_free(new_groups[j],S->groups[j].count*n,S->base->field);free(new_groups);}
    for(unsigned j=0;j<g*g;j++)b_theta4_vector_free(pairs[j],n,B->field);
    b_theta4_vector_free(raw_basis,g*n,B->field);b_theta4_vector_free(basis,g*n,B->field);
    b_theta4_vector_free(null,n,B->field);b_theta4_vector_free(raw,n,B->field);
    b_theta4_vector_free(new_generators,g*n,B->field);
    return ok;
}
