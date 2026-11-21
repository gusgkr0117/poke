#include "root-ec-tate.h"
static int rec(fmpz_t x,fmpz_t y,const a_ec_point *P,const a_ec_point B[2],unsigned bits,const fq_t A,const root_extension *E){
 if(bits==1){
  if(P->infinity){fmpz_zero(x);fmpz_zero(y);return 1;}
  for(int k=0;k<2;k++)if(a_ec_equal(P,B+k,E)){fmpz_set_ui(x,k==0);fmpz_set_ui(y,k==1);return 1;}
  a_ec_point S;a_ec_init(&S,E);a_ec_add(&S,B,B+1,A,E);int ok=a_ec_equal(P,&S,E);a_ec_clear(&S,E);if(ok){fmpz_one(x);fmpz_one(y);}return ok;
 }
 int ok=0;unsigned low=bits/2,high=bits-low;
 a_ec_point L[2],H[2],Q,R,T;for(int k=0;k<2;k++){a_ec_init(L+k,E);a_ec_init(H+k,E);}a_ec_init(&Q,E);a_ec_init(&R,E);a_ec_init(&T,E);
 fmpz_t scale,x0,y0,x1,y1;fmpz_init(scale);fmpz_init(x0);fmpz_init(y0);fmpz_init(x1);fmpz_init(y1);fmpz_one(scale);fmpz_mul_2exp(scale,scale,high);
 for(int k=0;k<2;k++)a_ec_mul(L+k,B+k,scale,A,E);a_ec_mul(&Q,P,scale,A,E);
 if(!rec(x0,y0,&Q,L,low,A,E))goto done;
 a_ec_mul(&Q,B,x0,A,E);a_ec_mul(&T,B+1,y0,A,E);a_ec_add(&Q,&Q,&T,A,E);a_ec_neg(&Q,&Q,E);a_ec_add(&R,P,&Q,A,E);
 fmpz_one(scale);fmpz_mul_2exp(scale,scale,low);for(int k=0;k<2;k++)a_ec_mul(H+k,B+k,scale,A,E);
 if(!rec(x1,y1,&R,H,high,A,E))goto done;
 fmpz_mul_2exp(x1,x1,low);fmpz_add(x,x0,x1);fmpz_mul_2exp(y1,y1,low);fmpz_add(y,y0,y1);ok=1;
done:
 fmpz_clear(y1);fmpz_clear(x1);fmpz_clear(y0);fmpz_clear(x0);fmpz_clear(scale);a_ec_clear(&T,E);a_ec_clear(&R,E);a_ec_clear(&Q,E);for(int k=0;k<2;k++){a_ec_clear(H+k,E);a_ec_clear(L+k,E);}return ok;
}
int root_ec_dlog2(fmpz_t x,fmpz_t y,const a_ec_point *P,const a_ec_point B[2],unsigned bits,const fq_t A,const root_extension *E){
 if(!bits||!a_ec_curve_valid(A,E)||!a_ec_on_curve(P,A,E)||!a_ec_on_curve(B,A,E)||!a_ec_on_curve(B+1,A,E))return 0;
 int ok=0;fmpz_t n,x0,y0;fmpz_init(n);fmpz_init(x0);fmpz_init(y0);fmpz_one(n);fmpz_mul_2exp(n,n,bits-1);
 a_ec_point top[2],Q,R;for(int k=0;k<2;k++)a_ec_init(top+k,E);a_ec_init(&Q,E);a_ec_init(&R,E);
 for(int k=0;k<2;k++){a_ec_mul(top+k,B+k,n,A,E);if(top[k].infinity)goto done;a_ec_add(&Q,top+k,top+k,A,E);if(!Q.infinity)goto done;}
 if(a_ec_equal(top,top+1,E))goto done;
 if(!rec(x0,y0,P,B,bits,A,E))goto done;
 a_ec_mul(&Q,B,x0,A,E);a_ec_mul(&R,B+1,y0,A,E);a_ec_add(&Q,&Q,&R,A,E);if(!a_ec_equal(&Q,P,E))goto done;
 fmpz_mul_2exp(n,n,1);if(fmpz_sgn(x0)<0||fmpz_sgn(y0)<0||fmpz_cmp(x0,n)>=0||fmpz_cmp(y0,n)>=0)goto done;
 fmpz_set(x,x0);fmpz_set(y,y0);ok=1;
done:
 a_ec_clear(&R,E);a_ec_clear(&Q,E);for(int k=0;k<2;k++)a_ec_clear(top+k,E);fmpz_clear(y0);fmpz_clear(x0);fmpz_clear(n);return ok;
}
int root_ec_maxorder_actions(fmpz_mat_struct out[4],const a_ec_point B[2],const a_ec_point half[2],unsigned bits,const root_extension *E){
 if(!bits)return 0;for(int k=0;k<4;k++)if(fmpz_mat_nrows(out+k)!=2||fmpz_mat_ncols(out+k)!=2)return 0;
 int ok=0;fmpz_mat_struct ans[4];for(int k=0;k<4;k++)fmpz_mat_init(ans+k,2,2);fmpz_mat_one(ans);
 fq_t A;fq_init(A,E->field);fq_zero(A,E->field);fmpz_t n,x,y;fmpz_init(n);fmpz_init(x);fmpz_init(y);fmpz_one(n);fmpz_mul_2exp(n,n,bits);
 fmpz c[4];for(int k=0;k<4;k++)fmpz_init(c+k);a_ec_point P;a_ec_init(&P,E);
 for(int k=1;k<4;k++){
  for(int j=0;j<4;j++)fmpz_set_ui(c+j,k==1?(j==1):k==2?(j==1||j==2):(j==0||j==3));
  for(int column=0;column<2;column++){
   if(a_ec_quaternion(&P,B+column,c,k==1?1:2,n,half+column,E)!=A_EC_OK||!root_ec_dlog2(x,y,&P,B,bits,A,E))goto done;
   fmpz_set(fmpz_mat_entry(ans+k,0,column),x);fmpz_set(fmpz_mat_entry(ans+k,1,column),y);
  }
 }
 for(int k=0;k<4;k++)fmpz_mat_set(out+k,ans+k);ok=1;
done:
 a_ec_clear(&P,E);for(int k=0;k<4;k++)fmpz_clear(c+k);fmpz_clear(y);fmpz_clear(x);fmpz_clear(n);fq_clear(A,E->field);for(int k=0;k<4;k++)fmpz_mat_clear(ans+k);return ok;
}
