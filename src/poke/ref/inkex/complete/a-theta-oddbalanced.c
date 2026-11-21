#include "a-theta-oddbalanced.h"
#include <stdlib.h>
#include <string.h>
void a_theta_oddbalanced_clear(a_theta_oddbalanced *S){if(!S)return;a_theta_oddfields_clear(&S->state);free(S->height);memset(S,0,sizeof(*S));}
int a_theta_oddbalanced_init(a_theta_oddbalanced *S,unsigned g,ulong ell,unsigned exponent,
 const root_extension *base,const root_extension *kernel_field,const fq_struct *null,
 const fq_struct *const *generators,const fq_struct *const *complements,const fq_struct *pairing4,
 const a_theta_payload_input *payloads,size_t groups){
 if(!S||S->ready||!exponent||!g||g>8||groups==SIZE_MAX||groups>SIZE_MAX/sizeof(a_theta_payload_input)-1)return 0;
 a_theta_payload_input *all=calloc(groups+1,sizeof(*all));if(!all)return 0;
 if(groups){if(!payloads){free(all);return 0;}memcpy(all,payloads,groups*sizeof(*all));}all[groups].field=kernel_field;all[groups].points=generators;all[groups].count=g;
 int good=a_theta_oddfields_init(&S->state,g,ell,exponent,base,kernel_field,null,generators,complements,pairing4,all,groups+1);free(all);if(!good)return 0;
 unsigned cap=1,t=exponent;while(t>1){t/=2;cap++;}S->height=calloc(cap,sizeof(*S->height));if(!S->height){a_theta_oddfields_clear(&S->state);return 0;}
 S->stack_group=groups;S->capacity=cap;S->depth=S->max_depth=1;S->height[0]=exponent;S->ready=1;return 1;
}
static int push(a_theta_oddbalanced *S){
 a_theta_oddfields *A=&S->state;a_theta_odd_payload_group *G=A->groups+S->stack_group;unsigned g=A->theta.g,h=S->height[S->depth-1],child=h/2,shift=h-child;size_t n=A->theta.count;
 if(h<2||S->depth>=S->capacity||G->count>SIZE_MAX/n-g)return 0;const root_extension *E=A->kernel_field;size_t old=G->count;
 fq_struct *next=b_theta4_vector_new((old+g)*n,E->field),*null=b_theta4_vector_new(n,E->field);b_theta4_ctx T={0};int ready=0,ok=0;if(!next||!null)goto done;
 for(size_t i=0;i<n;i++){if(E==A->base)fq_set(null+i,A->theta.null+i,E->field);else root_field_up(null+i,A->theta.null+i,&A->kernel_embedding);}
 if(!b_theta4_init(&T,g,null,E->field))goto done;ready=1;
 for(size_t i=0;i<old*n;i++)fq_set(next+i,G->points+i,E->field);
 for(unsigned j=0;j<g;j++){fq_struct *P=next+(old+j)*n;b_theta4_copy(P,G->points+(old-g+j)*n,&T);for(unsigned k=0;k<shift;k++)if(!b_theta4_scalar(P,P,A->ell,&T))goto done;}
 b_theta4_vector_free(G->points,old*n,E->field);G->points=next;next=NULL;G->count+=g;S->height[S->depth++]=child;if(S->depth>S->max_depth)S->max_depth=S->depth;S->frame_multiplications+=(unsigned long long)g*shift;ok=1;
done:if(ready)b_theta4_clear(&T);b_theta4_vector_free(null,n,E->field);b_theta4_vector_free(next,(old+g)*n,E->field);return ok;
}
int a_theta_oddbalanced_step(a_theta_oddbalanced *S){
 if(!S||!S->ready||!S->state.remaining||!S->depth)return 0;a_theta_oddfields *A=&S->state;unsigned g=A->theta.g;size_t n=A->theta.count;
 while(S->height[S->depth-1]>1)if(!push(S))return 0;
 a_theta_odd_payload_group *G=A->groups+S->stack_group;size_t after_count=G->count-g;fq_struct *after=after_count?b_theta4_vector_new(after_count*n,A->kernel_field->field):NULL;if(after_count&&!after)return 0;
 for(size_t i=0;i<g*n;i++)fq_set(A->generators+i,G->points+after_count*n+i,A->kernel_field->field);
 unsigned remaining=A->remaining;A->remaining=1;int ok=a_theta_oddfields_step(A);A->remaining=ok?remaining-1:remaining;
 if(!ok){b_theta4_vector_free(after,after_count*n,A->kernel_field->field);return 0;}G=A->groups+S->stack_group;
 for(size_t i=0;i<after_count*n;i++)fq_set(after+i,G->points+i,A->kernel_field->field);b_theta4_vector_free(G->points,G->count*n,A->kernel_field->field);G->points=after;G->count=after_count;
 S->frame_evaluations+=(unsigned long long)(S->depth+1)*g;for(unsigned i=0;i<S->depth;i++)S->height[i]--;S->depth--;return 1;
}
