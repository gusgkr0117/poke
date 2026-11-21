#ifndef INKE_RESEARCH_THETA_ODD_BASECHANGE_H
#define INKE_RESEARCH_THETA_ODD_BASECHANGE_H
#include "b-theta4-odd.h"
#include "root-field-embedding.h"
typedef struct {b_theta4_ctx source;b_theta4_odd_ctx quotient;int ready;} a_theta_odd_basechange;
/* Imports the already selected weighted chart and every affine lift literally.
 * No new root, kernel basis, complementary chart, or marking is selected. */
int a_theta_odd_basechange_init(a_theta_odd_basechange *out,const b_theta4_odd_ctx *Q,const root_field_embedding *L);
void a_theta_odd_basechange_clear(a_theta_odd_basechange *out);
#endif
