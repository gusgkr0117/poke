#ifndef INKE_RESEARCH_TATE_PROJECTOR_H
#define INKE_RESEARCH_TATE_PROJECTOR_H
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
typedef unsigned char root_tate_matrix[16][16];
/* Dimension-aware versions use the leading 2g-by-2g entries and zero
 * the unused storage. The two source blocks each contain g/2 elliptic
 * factors; g must be even, 2<=g<=8. Existing wrappers below use g=8. */
int root_tate_projector_dim(root_tate_matrix e,root_tate_matrix S,
 root_tate_matrix frame,const fmpz_mat_t P,const fmpz_mat_t eta,
 const fmpz_t n1,const fmpz_t N,unsigned g);
int root_tate_theta_frame_dim(root_tate_matrix out,const root_tate_matrix labels,
 const root_tate_matrix lattice_frame,unsigned g);
/* For a genuine elliptic graph, eta is its coherent 2x2 Tate matrix and
 * det(eta)=n1*n2 to the precision required below. P=[[n1,hat(eta)],
 * [eta,n2]], with scalar diagonal blocks; P is 4x4. */
int root_tate_projector_input_g2(fmpz_mat_t P,const fmpz_mat_t eta,
 const fmpz_t n1,const fmpz_t n2);
/* Source coordinates are four elliptic pairs in each of two blocks.  P is
 * the 16x16 numerator of the desired projector, eta is the 8x8 graph map.
 * N=n1+n2; n1 must be odd.  P and eta must come from ONE coherent splitting
 * of the endomorphism order modulo at least 2^(2*v2(N)+2).
 * Returns e=P/N and the target polarization S modulo4 in the explicit
 * lattice basis [[I/h,0],[G/h,I]], h=2^v2(N), G=eta/n1 modh.
 * Checks all divisions, alternating/unimodular S, self-adjoint idempotent e,
 * and constructs a symplectic frame with its first elliptic pair in im(e).
 * Output columns are alpha[0:8],beta[0:8]; e kills the other seven pairs.
 * All outputs remain unchanged on failure.  This does not perform a theta
 * transform or identify the current theta marking with this lattice basis. */
int root_tate_projector(root_tate_matrix e,root_tate_matrix S,
                         root_tate_matrix frame,const fmpz_mat_t P,
                         const fmpz_mat_t eta,const fmpz_t n1,const fmpz_t N);
/* Build the graph map and one projector numerator from integer four-square
 * matrices X,Y and a coherent 2x2 theta matrix.  d1*d2=det(theta),
 * XX^t=uI, YY^t=vI and XY=YX are caller promises checked by the caller. */
int root_tate_projector_input(fmpz_mat_t P,fmpz_mat_t eta,
                               const fmpz_mat_t X,const fmpz_mat_t Y,
                               const fmpz_mat_t theta,const fmpz_t d1,
                               const fmpz_t d2,unsigned factor);
/* labels gives current theta labels of the16 target-lattice basis4torsion
 * points. Convert an S-symplectic lattice frame to an actual Sp16(Z/4)
 * theta frame; normalize a possible common inverse fourth-root convention.
 * This requires all16 evaluated labels; abstract lattice data alone does
 * not determine the current theta marking. Outputs unchanged on failure. */
int root_tate_theta_frame(root_tate_matrix out,const root_tate_matrix labels,
                           const root_tate_matrix lattice_frame);
#endif
