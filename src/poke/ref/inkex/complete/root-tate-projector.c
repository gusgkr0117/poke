#include "root-tate-projector.h"
#include <string.h>
static int pairing(const unsigned char *a,const unsigned char *b,const root_tate_matrix S,int g){
 int n=2*g;
 int s=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)s+=a[i]*S[i][j]*b[j];return s&3;
}
static void mv(unsigned char *out,const root_tate_matrix A,const unsigned char *v,int g){
 int n=2*g;
 for(int i=0;i<n;i++){int z=0;for(int j=0;j<n;j++)z+=A[i][j]*v[j];out[i]=(unsigned char)(z&3);}
}
int root_tate_theta_frame_dim(root_tate_matrix out,const root_tate_matrix labels,const root_tate_matrix lattice_frame,unsigned dimension){
 if(dimension<2||dimension>8||(dimension&1))return 0;int g=(int)dimension,n=2*g;
 root_tate_matrix M={{0}};for(int i=0;i<n;i++)for(int j=0;j<n;j++){int s=0;for(int k=0;k<n;k++)s+=labels[i][k]*lattice_frame[k][j];M[i][j]=(unsigned char)(s&3);}
 int factor=0;for(int k=0;k<g;k++)factor+=M[k][0]*M[k+g][g]-M[k+g][0]*M[k][g];factor&=3;if(!(factor&1))return 0;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
  int s=0;for(int k=0;k<g;k++)s+=M[k][i]*M[k+g][j]-M[k+g][i]*M[k][j];
  int want=i<g&&j==i+g?factor:j<g&&i==j+g?-factor:0;if(((s-want)&3)!=0)return 0;
 }
 for(int i=0;i<n;i++)for(int j=g;j<n;j++)M[i][j]=(unsigned char)((factor*M[i][j])&3);memcpy(out,M,sizeof(M));return 1;
}
static int adapted_frame(root_tate_matrix F,const root_tate_matrix e,const root_tate_matrix S,int g){
 int n=2*g;
 unsigned char a[8][16]={{0}},b[8][16]={{0}},v[16][16],z[16];
 for(int k=0;k<g;k++){
  for(int j=0;j<n;j++){
   for(int i=0;i<n;i++)v[j][i]=(unsigned char)((k==0?e[i][j]:(i==j)-e[i][j])&3);
   for(int t=0;t<k;t++){
    int av=pairing(a[t],v[j],S,g),bv=pairing(b[t],v[j],S,g);
    for(int i=0;i<n;i++)v[j][i]=(unsigned char)((v[j][i]+bv*a[t][i]-av*b[t][i])&3);
   }
  }
  int x=-1,y=-1,p=0;for(int i=0;i<n&&x<0;i++)for(int j=i+1;j<n;j++){
   p=pairing(v[i],v[j],S,g);if(p&1){x=i;y=j;break;}
  }
  if(x<0)return 0;
  for(int i=0;i<n;i++){a[k][i]=v[x][i];b[k][i]=(unsigned char)((p*v[y][i])&3);}
 }
 memset(F,0,sizeof(root_tate_matrix));
 for(int k=0;k<g;k++)for(int i=0;i<n;i++){F[i][k]=a[k][i];F[i][k+g]=b[k][i];}
 /* Exact mod4 checks also certify unimodularity and the projector rank. */
 for(int i=0;i<g;i++)for(int j=0;j<g;j++){
  if(pairing(a[i],a[j],S,g)||pairing(b[i],b[j],S,g)||pairing(a[i],b[j],S,g)!=(i==j))return 0;
 }
 for(int k=0;k<g;k++){
  mv(z,e,a[k],g);for(int i=0;i<n;i++)if(z[i]!=(k==0?a[k][i]:0))return 0;
  mv(z,e,b[k],g);for(int i=0;i<n;i++)if(z[i]!=(k==0?b[k][i]:0))return 0;
 }
 return 1;
}
static int divide_mod4(unsigned char *out,const fmpz_t z,const fmpz_t divisor,ulong odd_inverse){
 if(!fmpz_divisible(z,divisor))return 0;fmpz_t q;fmpz_init(q);fmpz_divexact(q,z,divisor);
 *out=(unsigned char)((fmpz_fdiv_ui(q,4)*odd_inverse)&3);fmpz_clear(q);return 1;
}
int root_tate_projector_dim(root_tate_matrix out_e,root_tate_matrix out_S,root_tate_matrix out_F,
 const fmpz_mat_t P,const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t N,unsigned dimension){
 if(dimension<2||dimension>8||(dimension&1))return 0;int g=(int)dimension,n=2*g;
 if(fmpz_mat_nrows(P)!=n||fmpz_mat_ncols(P)!=n||fmpz_mat_nrows(eta)!=g||fmpz_mat_ncols(eta)!=g||!fmpz_is_odd(n1)||fmpz_sgn(N)<=0)return 0;
 int ok=0;root_tate_matrix e,S,F;memset(e,0,sizeof(e));memset(S,0,sizeof(S));
 fmpz_t h,h2,m,inv,t,z;fmpz_init(h);fmpz_init(h2);fmpz_init(m);fmpz_init(inv);fmpz_init(t);fmpz_init(z);
 ulong exponent=fmpz_val2(N);fmpz_one(h);fmpz_mul_2exp(h,h,exponent);fmpz_mul(h2,h,h);fmpz_divexact(m,N,h);ulong odd=fmpz_fdiv_ui(m,4);
 fmpz_mat_t G,A,B,C,D,BG,DG,GA,GBG,T;fmpz_mat_init(G,g,g);fmpz_mat_init(A,g,g);fmpz_mat_init(B,g,g);fmpz_mat_init(C,g,g);fmpz_mat_init(D,g,g);fmpz_mat_init(BG,g,g);fmpz_mat_init(DG,g,g);fmpz_mat_init(GA,g,g);fmpz_mat_init(GBG,g,g);fmpz_mat_init(T,g,g);
 if(exponent&&!fmpz_invmod(inv,n1,h))goto done;
 for(int i=0;i<g;i++)for(int j=0;j<g;j++){
  fmpz_mul(t,fmpz_mat_entry(eta,i,j),inv);fmpz_mod(fmpz_mat_entry(G,i,j),t,h);
  fmpz_set(fmpz_mat_entry(A,i,j),fmpz_mat_entry(P,i,j));fmpz_set(fmpz_mat_entry(B,i,j),fmpz_mat_entry(P,i,j+g));
  fmpz_set(fmpz_mat_entry(C,i,j),fmpz_mat_entry(P,i+g,j));fmpz_set(fmpz_mat_entry(D,i,j),fmpz_mat_entry(P,i+g,j+g));
 }
 fmpz_mat_mul(BG,B,G);fmpz_mat_mul(DG,D,G);fmpz_mat_mul(GA,G,A);fmpz_mat_mul(GBG,G,BG);fmpz_mat_mul(T,G,B);
 for(int i=0;i<g;i++)for(int j=0;j<g;j++){
  fmpz_add(z,fmpz_mat_entry(A,i,j),fmpz_mat_entry(BG,i,j));if(!divide_mod4(&e[i][j],z,h,odd))goto done;
  e[i][j+g]=(unsigned char)((fmpz_fdiv_ui(fmpz_mat_entry(B,i,j),4)*odd)&3);
  fmpz_add(z,fmpz_mat_entry(C,i,j),fmpz_mat_entry(DG,i,j));fmpz_sub(z,z,fmpz_mat_entry(GA,i,j));fmpz_sub(z,z,fmpz_mat_entry(GBG,i,j));if(!divide_mod4(&e[i+g][j],z,h2,odd))goto done;
  fmpz_sub(z,fmpz_mat_entry(D,i,j),fmpz_mat_entry(T,i,j));if(!divide_mod4(&e[i+g][j+g],z,h,odd))goto done;
  /* J is block diagonal with [[0,1],[-1,0]] on each elliptic factor. */
  fmpz_zero(z);for(int k=0;k<g/2;k++){
   fmpz_addmul(z,fmpz_mat_entry(G,2*k,i),fmpz_mat_entry(G,2*k+1,j));
   fmpz_submul(z,fmpz_mat_entry(G,2*k+1,i),fmpz_mat_entry(G,2*k,j));
  }
  int Jij=(i%2==0&&j==i+1)?1:(i%2==1&&j==i-1)?-1:0;fmpz_add_si(z,z,Jij);
  if(!divide_mod4(&S[i][j],z,h,odd))goto done; /* odd inverse=odd mod4 */
  fmpz_set(z,fmpz_mat_entry(G,j^1,i));if(!(j&1))fmpz_neg(z,z);S[i][j+g]=(unsigned char)((fmpz_fdiv_ui(z,4)*odd)&3);
  S[j+g][i]=(unsigned char)((-S[i][j+g])&3);
  fmpz_mul_si(z,N,Jij);S[i+g][j+g]=(unsigned char)fmpz_fdiv_ui(z,4);
 }
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){
  if((i==j&&S[i][j])||((S[i][j]+S[j][i])&3))goto done;
  int ee=0,se=0,ets=0;for(int k=0;k<n;k++){ee+=e[i][k]*e[k][j];se+=S[i][k]*e[k][j];ets+=e[k][i]*S[k][j];}
  if((ee&3)!=e[i][j]||((se-ets)&3))goto done;
 }
 if(!adapted_frame(F,e,S,g))goto done;
 memcpy(out_e,e,sizeof(e));memcpy(out_S,S,sizeof(S));memcpy(out_F,F,sizeof(F));ok=1;
done:
 fmpz_mat_clear(T);fmpz_mat_clear(GBG);fmpz_mat_clear(GA);fmpz_mat_clear(DG);fmpz_mat_clear(BG);fmpz_mat_clear(D);fmpz_mat_clear(C);fmpz_mat_clear(B);fmpz_mat_clear(A);fmpz_mat_clear(G);
 fmpz_clear(z);fmpz_clear(t);fmpz_clear(inv);fmpz_clear(m);fmpz_clear(h2);fmpz_clear(h);return ok;
}
int root_tate_projector(root_tate_matrix e,root_tate_matrix S,root_tate_matrix F,const fmpz_mat_t P,const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t N){return root_tate_projector_dim(e,S,F,P,eta,n1,N,8);}
int root_tate_theta_frame(root_tate_matrix out,const root_tate_matrix labels,const root_tate_matrix F){return root_tate_theta_frame_dim(out,labels,F,8);}
int root_tate_projector_input_g2(fmpz_mat_t P,const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t n2){
 if(fmpz_mat_nrows(P)!=4||fmpz_mat_ncols(P)!=4||fmpz_mat_nrows(eta)!=2||fmpz_mat_ncols(eta)!=2)return 0;
 fmpz_mat_zero(P);for(int i=0;i<2;i++){fmpz_set(fmpz_mat_entry(P,i,i),n1);fmpz_set(fmpz_mat_entry(P,i+2,i+2),n2);}
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){fmpz_set(fmpz_mat_entry(P,i+2,j),fmpz_mat_entry(eta,i,j));fmpz_set(fmpz_mat_entry(P,i,j+2),fmpz_mat_entry(eta,1-j,1-i));if(i!=j)fmpz_neg(fmpz_mat_entry(P,i,j+2),fmpz_mat_entry(P,i,j+2));}return 1;
}
int root_tate_projector_input(fmpz_mat_t P,fmpz_mat_t eta,const fmpz_mat_t X,const fmpz_mat_t Y,const fmpz_mat_t theta,const fmpz_t d1,const fmpz_t d2,unsigned r){
 if(r>=4||fmpz_mat_nrows(P)!=16||fmpz_mat_ncols(P)!=16||fmpz_mat_nrows(eta)!=8||fmpz_mat_ncols(eta)!=8||fmpz_mat_nrows(X)!=4||fmpz_mat_ncols(X)!=4||fmpz_mat_nrows(Y)!=4||fmpz_mat_ncols(Y)!=4||fmpz_mat_nrows(theta)!=2||fmpz_mat_ncols(theta)!=2)return 0;
 fmpz_t t;fmpz_init(t);fmpz_mat_t YX,hat;fmpz_mat_init(YX,4,4);fmpz_mat_init(hat,2,2);fmpz_mat_mul(YX,Y,X);
 fmpz_set(fmpz_mat_entry(hat,0,0),fmpz_mat_entry(theta,1,1));fmpz_set(fmpz_mat_entry(hat,1,1),fmpz_mat_entry(theta,0,0));fmpz_neg(fmpz_mat_entry(hat,0,1),fmpz_mat_entry(theta,0,1));fmpz_neg(fmpz_mat_entry(hat,1,0),fmpz_mat_entry(theta,1,0));fmpz_mat_zero(P);
 for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int a=0;a<2;a++)for(int b=0;b<2;b++){
  fmpz_mul(fmpz_mat_entry(eta,2*i+a,2*j+b),fmpz_mat_entry(YX,i,j),fmpz_mat_entry(theta,a,b));
  if(a==b){
   fmpz_mul(t,fmpz_mat_entry(X,r,i),fmpz_mat_entry(X,r,j));fmpz_mul(fmpz_mat_entry(P,2*i+a,2*j+b),t,d1);
   fmpz_mul(t,fmpz_mat_entry(Y,i,r),fmpz_mat_entry(Y,j,r));fmpz_mul(fmpz_mat_entry(P,8+2*i+a,8+2*j+b),t,d2);
  }
  fmpz_mul(t,fmpz_mat_entry(X,r,i),fmpz_mat_entry(Y,j,r));fmpz_mul(fmpz_mat_entry(P,2*i+a,8+2*j+b),t,fmpz_mat_entry(hat,a,b));
  fmpz_mul(t,fmpz_mat_entry(Y,i,r),fmpz_mat_entry(X,r,j));fmpz_mul(fmpz_mat_entry(P,8+2*i+a,2*j+b),t,fmpz_mat_entry(theta,a,b));
 }
 fmpz_mat_clear(hat);fmpz_mat_clear(YX);fmpz_clear(t);return 1;
}
