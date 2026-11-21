#include "root-quaternion-splitting.h"
static void reduce(fmpz_mat_t A,const fmpz_t m){for(slong i=0;i<fmpz_mat_nrows(A);i++)for(slong j=0;j<fmpz_mat_ncols(A);j++)fmpz_mod(fmpz_mat_entry(A,i,j),fmpz_mat_entry(A,i,j),m);}
static int inverse_mod(fmpz_mat_t inv,const fmpz_mat_t A,const fmpz_t m){
 fmpz_t den,x;fmpz_init(den);fmpz_init(x);int ok=fmpz_mat_inv(inv,den,A)&&fmpz_invmod(x,den,m);
 if(ok){fmpz_mat_scalar_mul_fmpz(inv,inv,x);reduce(inv,m);}fmpz_clear(x);fmpz_clear(den);return ok;
}
int root_quaternion_splitting(fmpz_mat_struct out[4],const fmpz_t p,unsigned bits){
 if(bits<1||fmpz_fdiv_ui(p,8)!=7)return 0;for(int i=0;i<4;i++)if(fmpz_mat_nrows(out+i)!=2||fmpz_mat_ncols(out+i)!=2)return 0;
 fmpz_t s,z,t,m;fmpz_init(s);fmpz_init(z);fmpz_init(t);fmpz_init(m);fmpz_one(s);fmpz_one(m);fmpz_mul_2exp(m,m,bits);
 /* Half-integral basis elements have norm denominator4.  The root equation
  * therefore needs TWO guard bits, not only the one used for entry division. */
 for(unsigned k=3;k<bits+2;k++){fmpz_mul(t,s,s);fmpz_add(t,t,p);if(fmpz_tstbit(t,k)){fmpz_one(z);fmpz_mul_2exp(z,z,k-1);fmpz_add(s,s,z);}}
 fmpz_mat_t j;fmpz_mat_init(j,2,2);fmpz_set(fmpz_mat_entry(j,0,0),s);fmpz_neg(fmpz_mat_entry(j,1,0),s);fmpz_neg(fmpz_mat_entry(j,1,1),s);
 fmpz_mat_one(out);fmpz_set_si(fmpz_mat_entry(out+1,0,0),-1);fmpz_set_si(fmpz_mat_entry(out+1,0,1),-2);fmpz_one(fmpz_mat_entry(out+1,1,0));fmpz_one(fmpz_mat_entry(out+1,1,1));
 fmpz_mat_add(out+2,out+1,j);fmpz_mat_scalar_divexact_ui(out+2,out+2,2);
 fmpz_mat_mul(out+3,out+1,j);fmpz_mat_add(out+3,out+3,out);fmpz_mat_scalar_divexact_ui(out+3,out+3,2);
 for(int i=0;i<4;i++)reduce(out+i,m);fmpz_mat_clear(j);fmpz_clear(m);fmpz_clear(t);fmpz_clear(z);fmpz_clear(s);return 1;
}
int root_quaternion_splitting_align(fmpz_mat_struct out[4],const fmpz_mat_struct known[4],const fmpz_t p,unsigned known_bits,unsigned high_bits){
 if(known_bits<1||high_bits<known_bits)return 0;for(int i=0;i<4;i++)if(fmpz_mat_nrows(out+i)!=2||fmpz_mat_ncols(out+i)!=2||fmpz_mat_nrows(known+i)!=2||fmpz_mat_ncols(known+i)!=2)return 0;
 int ok=0;fmpz_t m,q,z;fmpz_init(m);fmpz_init(q);fmpz_init(z);fmpz_one(m);fmpz_mul_2exp(m,m,high_bits);fmpz_one(q);fmpz_mul_2exp(q,q,known_bits);
 fmpz_mat_struct standard[4],answer[4];for(int i=0;i<4;i++){fmpz_mat_init(standard+i,2,2);fmpz_mat_init(answer+i,2,2);}
 fmpz_mat_t V,Vinv,E11,E21,U,Uinv,T;fmpz_mat_init(V,4,4);fmpz_mat_init(Vinv,4,4);fmpz_mat_init(E11,2,2);fmpz_mat_init(E21,2,2);fmpz_mat_init(U,2,2);fmpz_mat_init(Uinv,2,2);fmpz_mat_init(T,2,2);
 if(!root_quaternion_splitting(standard,p,high_bits))goto done;
 for(int k=0;k<4;k++)for(int a=0;a<2;a++)for(int b=0;b<2;b++)fmpz_mod(fmpz_mat_entry(V,2*a+b,k),fmpz_mat_entry(standard+k,a,b),q);
 if(!inverse_mod(Vinv,V,q))goto done;
 for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int k=0;k<4;k++){
  fmpz_addmul(fmpz_mat_entry(E11,a,b),fmpz_mat_entry(Vinv,k,0),fmpz_mat_entry(known+k,a,b));
  fmpz_addmul(fmpz_mat_entry(E21,a,b),fmpz_mat_entry(Vinv,k,2),fmpz_mat_entry(known+k,a,b));
 }
 reduce(E11,q);reduce(E21,q);int col=-1;
 for(int b=0;b<2&&col<0;b++)for(int a=0;a<2;a++)if(fmpz_is_odd(fmpz_mat_entry(E11,a,b))){col=b;break;}
 if(col<0)goto done;
 for(int a=0;a<2;a++)fmpz_set(fmpz_mat_entry(U,a,0),fmpz_mat_entry(E11,a,col));
 for(int a=0;a<2;a++)for(int b=0;b<2;b++)fmpz_addmul(fmpz_mat_entry(U,a,1),fmpz_mat_entry(E21,a,b),fmpz_mat_entry(U,b,0));
 reduce(U,q);if(!inverse_mod(Uinv,U,m))goto done;
 for(int k=0;k<4;k++){
  fmpz_mat_mul(T,U,standard+k);fmpz_mat_mul(answer+k,T,Uinv);reduce(answer+k,m);
  for(int a=0;a<2;a++)for(int b=0;b<2;b++){
   fmpz_sub(z,fmpz_mat_entry(answer+k,a,b),fmpz_mat_entry(known+k,a,b));if(!fmpz_divisible(z,q))goto done;
  }
 }
 for(int k=0;k<4;k++)fmpz_mat_set(out+k,answer+k);ok=1;
done:
 fmpz_mat_clear(T);fmpz_mat_clear(Uinv);fmpz_mat_clear(U);fmpz_mat_clear(E21);fmpz_mat_clear(E11);fmpz_mat_clear(Vinv);fmpz_mat_clear(V);for(int i=0;i<4;i++){fmpz_mat_clear(answer+i);fmpz_mat_clear(standard+i);}fmpz_clear(z);fmpz_clear(q);fmpz_clear(m);return ok;
}
int root_quaternion_splitting_element(fmpz_mat_t out,const fmpz *coeff,unsigned den,const fmpz_mat_struct basis[4],unsigned bits){
 if((den!=1&&den!=2)||bits<1||fmpz_mat_nrows(out)!=2||fmpz_mat_ncols(out)!=2)return 0;
 fmpz c[4];for(int i=0;i<4;i++)fmpz_init(c+i);fmpz_sub(c,coeff,coeff+3);fmpz_sub(c+1,coeff+1,coeff+2);fmpz_mul_2exp(c+2,coeff+2,1);fmpz_mul_2exp(c+3,coeff+3,1);
 int ok=1;for(int k=0;k<4;k++){if(den==2&&fmpz_is_odd(c+k)){ok=0;break;}fmpz_divexact_ui(c+k,c+k,den);}
 if(ok){fmpz_t m;fmpz_init(m);fmpz_one(m);fmpz_mul_2exp(m,m,bits);fmpz_mat_t A;fmpz_mat_init(A,2,2);
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int k=0;k<4;k++)fmpz_addmul(fmpz_mat_entry(A,a,b),c+k,fmpz_mat_entry(basis+k,a,b));reduce(A,m);fmpz_mat_set(out,A);fmpz_mat_clear(A);fmpz_clear(m);
 }
 for(int i=0;i<4;i++)fmpz_clear(c+i);return ok;
}
