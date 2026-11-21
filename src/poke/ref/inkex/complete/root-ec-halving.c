#include "root-ec-halving.h"
/* Bekker--Zarhin, arXiv:1702.02255, Section2 formula(6), specialized to
 * roots0,i,-i. Choose r0*r1*r2=-y; put s1=sum r_i,s2=sum r_i*r_j.
 * Q=(x+s2,-y-s1*s2). The signed doubling check is unconditional. */
int root_ec_half_e0(a_ec_point *out,const a_ec_point *P,const root_extension *E){
 fq_t A,r[3],s1,s2,t,product;fq_init(A,E->field);fq_zero(A,E->field);for(int k=0;k<3;k++)fq_init(r[k],E->field);fq_init(s1,E->field);fq_init(s2,E->field);fq_init(t,E->field);fq_init(product,E->field);a_ec_point Q,twice;a_ec_init(&Q,E);a_ec_init(&twice,E);int ok=0;
 if(!a_ec_on_curve(P,A,E))goto done;if(P->infinity){a_ec_infinity(&Q,E);ok=1;goto success;}
 if(!root_extension_sqrt(r[0],P->x,E))goto done;fq_sub(t,P->x,E->i,E->field);if(!root_extension_sqrt(r[1],t,E))goto done;fq_add(t,P->x,E->i,E->field);if(!root_extension_sqrt(r[2],t,E))goto done;
 fq_mul(product,r[0],r[1],E->field);fq_mul(product,product,r[2],E->field);fq_neg(t,P->y,E->field);if(!fq_equal(product,t,E->field)){fq_neg(r[0],r[0],E->field);fq_neg(product,product,E->field);}if(!fq_equal(product,t,E->field))goto done;
 fq_add(s1,r[0],r[1],E->field);fq_add(s1,s1,r[2],E->field);fq_mul(s2,r[0],r[1],E->field);fq_mul(t,r[0],r[2],E->field);fq_add(s2,s2,t,E->field);fq_mul(t,r[1],r[2],E->field);fq_add(s2,s2,t,E->field);
 fq_add(Q.x,P->x,s2,E->field);fq_mul(Q.y,s1,s2,E->field);fq_add(Q.y,Q.y,P->y,E->field);fq_neg(Q.y,Q.y,E->field);Q.infinity=0;
 if(!a_ec_on_curve(&Q,A,E))goto done;a_ec_add(&twice,&Q,&Q,A,E);if(!a_ec_equal(&twice,P,E))goto done;ok=1;
success:if(ok)a_ec_set(out,&Q,E);
done:a_ec_clear(&twice,E);a_ec_clear(&Q,E);fq_clear(product,E->field);fq_clear(t,E->field);fq_clear(s2,E->field);fq_clear(s1,E->field);for(int k=0;k<3;k++)fq_clear(r[k],E->field);fq_clear(A,E->field);return ok;
}
int root_ec_refine_actions(fmpz_mat_struct out[4],const fmpz_mat_struct known[4],
 const a_ec_point B[2],const a_ec_point half[2],unsigned low,unsigned bits,const root_extension *E){
 if(!low||bits<=low)return 0;for(int k=0;k<4;k++)if(fmpz_mat_nrows(out+k)!=2||fmpz_mat_ncols(out+k)!=2||fmpz_mat_nrows(known+k)!=2||fmpz_mat_ncols(known+k)!=2)return 0;
 int ok=0;fmpz_mat_struct ans[4];for(int k=0;k<4;k++)fmpz_mat_init(ans+k,2,2);fmpz_mat_one(ans);fq_t A;fq_init(A,E->field);fq_zero(A,E->field);
 fmpz_t n,mod,x,y,x0,y0;fmpz_init(n);fmpz_init(mod);fmpz_init(x);fmpz_init(y);fmpz_init(x0);fmpz_init(y0);fmpz_one(n);fmpz_mul_2exp(n,n,bits);fmpz_one(mod);fmpz_mul_2exp(mod,mod,low);
 a_ec_point H[2],P,R,T;for(int k=0;k<2;k++){a_ec_init(H+k,E);a_ec_mul(H+k,B+k,mod,A,E);}a_ec_init(&P,E);a_ec_init(&R,E);a_ec_init(&T,E);fmpz c[4];for(int k=0;k<4;k++)fmpz_init(c+k);
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){fmpz_mod(x0,fmpz_mat_entry(known,i,j),mod);if(!fmpz_equal_ui(x0,i==j))goto done;}
 for(int k=1;k<4;k++){
  for(int j=0;j<4;j++)fmpz_set_ui(c+j,k==1?(j==1):k==2?(j==1||j==2):(j==0||j==3));
  for(int col=0;col<2;col++){
   if(a_ec_quaternion(&P,B+col,c,k==1?1:2,n,half+col,E)!=A_EC_OK)goto done;
   fmpz_mod(x0,fmpz_mat_entry(known+k,0,col),mod);fmpz_mod(y0,fmpz_mat_entry(known+k,1,col),mod);
   a_ec_mul(&R,B,x0,A,E);a_ec_mul(&T,B+1,y0,A,E);a_ec_add(&R,&R,&T,A,E);a_ec_neg(&R,&R,E);a_ec_add(&R,&P,&R,A,E);
   if(!root_ec_dlog2(x,y,&R,H,bits-low,A,E))goto done;
   fmpz_mul_2exp(x,x,low);fmpz_add(fmpz_mat_entry(ans+k,0,col),x0,x);fmpz_mul_2exp(y,y,low);fmpz_add(fmpz_mat_entry(ans+k,1,col),y0,y);
  }
 }
 for(int k=0;k<4;k++)fmpz_mat_set(out+k,ans+k);ok=1;
done:for(int k=0;k<4;k++)fmpz_clear(c+k);a_ec_clear(&T,E);a_ec_clear(&R,E);a_ec_clear(&P,E);for(int k=0;k<2;k++)a_ec_clear(H+k,E);fmpz_clear(y0);fmpz_clear(x0);fmpz_clear(y);fmpz_clear(x);fmpz_clear(mod);fmpz_clear(n);fq_clear(A,E->field);for(int k=0;k<4;k++)fmpz_mat_clear(ans+k);return ok;
}
