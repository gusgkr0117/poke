#!/usr/bin/env python3
"""Offline, deterministic INKE manuscript constants (standard-library only).

The resulting scheme is entirely C.  This program is not called by Keygen.
Half-integral endomorphisms are evaluated on a half of the INPUT point over
Fp4; halving their output instead would choose an unrelated 2-torsion lift.
Every action matrix is checked against points and the order multiplication.
"""
import argparse
import functools
import hashlib
import json
import math
import shutil
from pathlib import Path
import time

PARAMS = {1: (256, 259, 162, 162, 31), 3: (384, 386, 243, 244, 23),
          5: (512, 517, 324, 327, 1)}


class F2(tuple):
    p = 0
    def __new__(cls, a=0, b=0):
        if isinstance(a, F2): return a
        return tuple.__new__(cls, (a % cls.p, b % cls.p))
    def __add__(a, b):
        b = F2(b); return F2(a[0]+b[0], a[1]+b[1])
    __radd__ = __add__
    def __neg__(a): return F2(-a[0], -a[1])
    def __sub__(a, b): return a + -F2(b)
    def __rsub__(a, b): return F2(b) + -a
    def __mul__(a, b):
        b = F2(b); return F2(a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    __rmul__ = __mul__
    def inv(a):
        n = pow((a[0]*a[0]+a[1]*a[1]) % a.p, -1, a.p)
        return F2(a[0]*n, -a[1]*n)
    def __truediv__(a, b): return a * F2(b).inv()
    def __pow__(a, n):
        out = F2(1)
        while n:
            if n & 1: out *= a
            a *= a; n >>= 1
        return out
    def conj(a): return F2(a[0], -a[1])
    def sqrt(a):
        if a == F2(0): return F2(0)
        if a[1] == 0:
            z = pow(a[0], (a.p+1)//4, a.p)
            if z*z % a.p == a[0]: return F2(z)
            return F2(0, pow(-a[0] % a.p, (a.p+1)//4, a.p))
        norm = (a[0]*a[0]+a[1]*a[1]) % a.p
        n = pow(norm, (a.p+1)//4, a.p)
        if n*n % a.p != norm: return None
        for s in (n, -n):
            r = (a[0]+s) * ((a.p+1)//2) % a.p
            z = pow(r, (a.p+1)//4, a.p)
            if z and z*z % a.p == r:
                root = F2(z, a[1]*pow(2*z,-1,a.p))
                assert root*root == a
                return root
        raise AssertionError('Fp2 sqrt')


class F4(tuple):
    delta = None
    frob_v = None
    def __new__(cls, a=0, b=0):
        if isinstance(a,F4): return a
        return tuple.__new__(cls, (F2(a),F2(b)))
    def __add__(a,b):
        b=F4(b); return F4(a[0]+b[0],a[1]+b[1])
    __radd__=__add__
    def __neg__(a): return F4(-a[0],-a[1])
    def __sub__(a,b): return a+-F4(b)
    def __rsub__(a,b): return F4(b)+-a
    def __mul__(a,b):
        b=F4(b); return F4(a[0]*b[0]+a.delta*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    __rmul__=__mul__
    def inv(a):
        n=(a[0]*a[0]-a.delta*a[1]*a[1]).inv()
        return F4(a[0]*n,-a[1]*n)
    def __truediv__(a,b): return a*F4(b).inv()
    def frobenius(a): return F4(a[0].conj(),a[1].conj()*a.frob_v)
    @staticmethod
    def sqrt2(a):
        r=a.sqrt()
        if r is not None: return F4(r)
        r=(a/F4.delta).sqrt()
        assert r is not None
        return F4(0,r)


def neg(P): return None if P is None else (P[0],-P[1])
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x == u:
        if y != v or y == type(y)(0): return None
        s=(3*x*x+1)/(2*y)
    else: s=(v-y)/(u-x)
    z=s*s-x-u
    return z,s*(x-z)-y
def mul(P,n):
    if n < 0: return mul(neg(P),-n)
    out=None
    while n:
        if n&1: out=add(out,P)
        P=add(P,P); n>>=1
    return out
def ispoint(P): return P is None or P[1]*P[1]==P[0]*P[0]*P[0]+P[0]
def aut_i(P):
    if P is None:return P
    return -P[0], type(P[0])(F2(0,1))*P[1]
def frob(P):
    if P is None:return P
    return tuple(x.conj() if isinstance(x,F2) else x.frobenius() for x in P)
def half(P):
    x,y=P; ii=F2(0,1)
    r=[F4.sqrt2(x-e) for e in (F2(0),ii,-ii)]
    if r[0]*r[1]*r[2] != F4(-y): r[0]=-r[0]
    assert r[0]*r[1]*r[2] == F4(-y)
    s=r[0]*r[1]+r[1]*r[2]+r[2]*r[0]
    H=F4(x)+s,-F4(y)-(r[0]+r[1]+r[2])*s
    assert ispoint(H) and add(H,H)==(F4(x),F4(y))
    return H
def half_i_j(P):
    H=half(P)
    out=add(aut_i(H),frob(H))
    if out is None:return None
    assert out[0][1]==out[1][1]==F2(0)
    out=out[0][0],out[1][0]
    assert ispoint(out) and mul(out,2)==add(aut_i(P),frob(P))
    # Any input half yields the same result: i+j kills the entire 2-torsion.
    for x in (F2(0),F2(0,1),F2(0,-1)):
        H2=add(H,(F4(x),F4(0)))
        assert add(aut_i(H2),frob(H2))==(F4(out[0]),F4(out[1]))
    return out


def basis(ell,e,p):
    N=ell**e; top=N//ell
    for counter in range(1,10000):
        x=F2(counter,1); y=(x*x*x+x).sqrt()
        if y is None:continue
        P=mul((x,y),(p+1)//N); Q=frob(P)
        A=mul(P,top); B=mul(Q,top)
        if A is None or B is None:continue
        if B in [mul(A,k) for k in range(ell)]:continue
        assert mul(P,N) is None and mul(Q,N) is None
        print(f'  basis {ell}^{e}: x counter {counter}',flush=True)
        return P,Q
    raise RuntimeError('basis search')


def dlog(P,Q,R,ell,e):
    Ps=[P];Qs=[Q]
    for _ in range(e):Ps.append(mul(Ps[-1],ell));Qs.append(mul(Qs[-1],ell))
    table={add(mul(Ps[e-1],a),mul(Qs[e-1],b)):(a,b)
           for a in range(ell) for b in range(ell)}
    assert len(table)==ell*ell
    def rec(R,start,length):
        if length==1:return table[R]
        lo=length//2;hi=length-lo
        a,b=rec(mul(R,ell**hi),start+hi,lo)
        rest=add(add(R,neg(mul(Ps[start],a))),neg(mul(Qs[start],b)))
        c,d=rec(rest,start+lo,hi)
        return a+ell**lo*c,b+ell**lo*d
    a,b=rec(R,0,e)
    assert add(mul(P,a),mul(Q,b))==R
    return a,b


def matadd(A,B,n):return [[(A[i][j]+B[i][j])%n for j in range(2)]for i in range(2)]
def matmul(A,B,n):return [[sum(A[i][k]*B[k][j]for k in range(2))%n for j in range(2)]for i in range(2)]
def matscale(A,s,n):return [[x*s%n for x in row]for row in A]
def action(P,Q,fun,ell,e):
    col=[dlog(P,Q,fun(T),ell,e)for T in(P,Q)]
    return [[col[j][i]for j in range(2)]for i in range(2)]
def check_actions(I,J,H,n,p):
    U=[[1,0],[0,1]];Z=[[0,0],[0,0]]
    assert matmul(I,I,n)==matscale(U,-1,n)
    assert matmul(J,J,n)==matscale(U,-p,n)
    assert matadd(matmul(I,J,n),matmul(J,I,n),n)==Z
    assert matscale(H,2,n)==matadd(I,J,n)
    assert matmul(H,H,n)==matscale(U,-(p+1)//4,n)
    assert matadd(matmul(I,H,n),matmul(H,I,n),n)==matscale(U,-1,n)
    G=matadd(U,matmul(I,H,n),n)
    assert matmul(G,G,n)==matadd(G,matscale(U,-(p+1)//4,n),n)
    for A,norm,trace in((I,1,0),(J,p,0),(H,(p+1)//4,0),(G,(p+1)//4,1)):
        assert (A[0][0]*A[1][1]-A[0][1]*A[1][0])%n==norm%n
        assert (A[0][0]+A[1][1])%n==trace%n
    return G


def run(level,out):
    a,f2,b,f3,c=PARAMS[level]; p=2**f2*3**f3*c-1
    if f2 < a + 2 or not 4*3**(2*b) < p < 3*2**(2*a-4)*3**b:
        raise ValueError('The field must support the auxiliary lifts and both manuscript bounds')
    F2.p=p
    for k in range(1,100):
        d=F2(k,1)
        if d.sqrt() is None:break
    F4.delta=d;F4.frob_v=d**((p-1)//2)
    assert F4.frob_v*F4.frob_v==d.conj()/d
    print(f'level {level}: {p.bit_length()}-bit p, Fp4 nonresidue {k}+i',flush=True)
    data={'level':level,'a':a,'f2':f2,'b':b,'f3':f3,'cofactor':c,'p':p,'basis':{},'actions':{}}
    acts=[]
    for ell,e in((2,f2),(3,b)):
        P,Q=basis(ell,e,p); n=ell**e
        data['basis'][str(ell)]=[P,Q,add(P,neg(Q))]
        print(f'  action i on {ell}-power torsion',flush=True)
        I=action(P,Q,aut_i,ell,e);J=[[0,1],[1,0]]
        print(f'  action (i+j)/2 on {ell}-power torsion',flush=True)
        if ell==2:H=action(P,Q,half_i_j,ell,e)
        else:H=matscale(matadd(I,J,n),pow(2,-1,n),n)
        G=check_actions(I,J,H,n,p)
        # Point-level checks for the other half-integral basis element.
        for col,T in enumerate((P,Q)):
            actual=add(T,aut_i(half_i_j(T))) if ell==2 else mul(add(T,aut_i(frob(T))),pow(2,-1,n))
            assert add(mul(P,G[0][col]),mul(Q,G[1][col]))==actual
        acts.append((n,[I,J,matmul(I,J,n),I,H,G]))
        print(f'  all point actions and maximal-order relations verified modulo {ell}^{e}',flush=True)
    n,m=acts[0][0],acts[1][0]
    for index,name in enumerate(('I','J','K','GEN2','GEN3','GEN4')):
        A=acts[0][1][index];B=acts[1][1][index]
        data['actions'][name]=[[A[i][j]+n*((B[i][j]-A[i][j])*pow(n,-1,m)%m)for j in range(2)]for i in range(2)]
    check_actions(data['actions']['I'],data['actions']['J'],data['actions']['GEN3'],n*m,p)
    out.mkdir(parents=True,exist_ok=True)
    (out/'constants.json').write_text(json.dumps(data,indent=2)+'\n')
    print(f'  wrote validated {out}/constants.json',flush=True)
    return data


def emit(data,out):
    """Emit read-only C constants; no executable generator is linked."""
    level=data['level'];a=data['a'];f=data['f2'];b=data['b'];p=data['p']
    F2.p=p; w=(p.bit_length()+63)//64; R=pow(2,64*w,p)
    inc=out/'include';inc.mkdir(exist_ok=True)
    def put(path,text):path.write_text(text+'\n',encoding='ascii')
    def limbs(x,n=w):return '{'+','.join(f'UINT64_C(0x{(x>>(64*i))&((1<<64)-1):016x})'for i in range(n))+'}'
    def ibz(x):
        n=(abs(x).bit_length()+63)//64
        return '{{._mp_alloc=0,._mp_size='+str(-n if x<0 else n)+',._mp_d=(mp_limb_t[])'+limbs(abs(x),max(n,1))+'}}'
    def imat(M):return '{'+','.join('{'+','.join(ibz(v)for v in row)+'}'for row in M)+'}'
    def fe(z):return '{'+limbs(z[0]*R%p)+','+limbs(z[1]*R%p)+'}'
    def point(P):return '{'+fe(P[0])+','+fe((1,0))+'}'
    def cbasis(B):return '{'+','.join(point(P)for P in B)+'}'
    banner='/* Generated by scripts/inke_precompute.py. Do not edit. */\n'
    put(inc/'fp_constants.h',banner+f'#define NWORDS_FIELD {w}\n#define NWORDS_ORDER {w}\n#define BITS {p.bit_length()}\n#define LOG2P {(p.bit_length()-1).bit_length()}\n#include <inke_field_constants.h>')
    put(inc/'inke_field_constants.h',banner+f'#define INKE_FIELD_MODULUS {limbs(p)}\n#define INKE_FIELD_R {limbs(R)}\n#define INKE_FIELD_R2 {limbs(R*R%p)}\n#define INKE_PUBLIC_TWO_EXPONENT {a}\n#define INKE_MANUSCRIPT_PARAMETERS 1')

    objects={}
    for name,which in(('EVEN','2'),('THREE','3'),('ODD_PLUS','3'),('COMMITMENT_PLUS','3'),('CHALLENGE','2')):
        objects['BASIS_'+name]=('ec_basis_t',cbasis(data['basis'][which]))
    # The twist subgroup has order one and is unused by INKE.
    identity_point='{'+fe((1,0))+','+fe((0,0))+'}'
    for name in('BASIS_ODD_MINUS','BASIS_COMMITMENT_MINUS'):
        objects[name]=('ec_basis_t','{'+','.join([identity_point]*3)+'}')
    objects['CURVE_E0']=('ec_curve_t','{.A='+fe((0,0))+',.C='+fe((1,0))+'}')
    objects['CURVE_E0_A24']=('ec_point_t','{'+fe((2,0))+','+fe((4,0))+'}')
    for name,M in data['actions'].items():objects['ACTION_'+name]=('ibz_mat_2x2_t',imat(M))
    # A primitive generator for the ideal annihilating the first 3^b-basis point.
    n=3**b;I=data['actions']['I'];x=-I[0][0]%n;z=-I[1][0]%n
    for attempt in range(12):
        xx=x+(attempt%3)*n;zz=z+(attempt//3)*n
        if math.gcd(xx*xx+1+p*zz*zz,n*n)==n:
            x,z=xx,zz;break
    else:raise AssertionError('primitive commitment ideal generator')
    objects['COMMITMENT_IDEAL_UNDISTORTED_GEN']=('quat_alg_elem_t','{'+ibz(1)+',{'+','.join(ibz(v)for v in(x,1,z,0))+'}}')
    objects['COMMITMENT_IDEAL_DISTORTION_ENDO']=('quat_alg_elem_t','{'+ibz(1)+',{'+','.join(ibz(v)for v in(0,0,1,0))+'}}')
    put(inc/'endomorphism_action.h',banner+'#include <intbig.h>\n#include <ec.h>\n#include <quaternion.h>\n'+'\n'.join(f'extern const {ty} {name};'for name,(ty,val)in objects.items()))
    put(out/'endomorphism_action.c',banner+'#include <stdint.h>\n#include <endomorphism_action.h>\n#if GMP_LIMB_BITS != 64\n#error "Generated constants require 64-bit GMP limbs"\n#endif\n'+'\n'.join(f'const {ty} {name}={val};'for name,(ty,val)in objects.items()))

    integer_values={'CHARACTERISTIC':p,'TORSION_ODD':n,'TORSION_ODD_PLUS':n,'TORSION_ODD_MINUS':1,
                    'TORSION_PLUS_2POWER':2**f,'TORSION_PLUS_3POWER':n,'TORSION_PLUS_23POWER':2**f*n,
                    'DEGREE_COMMITMENT':n,'DEGREE_COMMITMENT_PLUS':n,'DEGREE_COMMITMENT_MINUS':1,'DEGREE_CHALLENGE':2**f}
    h=banner+'#include <intbig.h>\n#include <stddef.h>\n'
    for name,value in(('2POWER',2**f),('3POWER',n),('23POWER',2**f*n)):
        h+=f'#define TORSION_{name}_BYTES {(value.bit_length()+7)//8}\n'
    c=banner+'#include <stdint.h>\n#include <torsion_constants.h>\n'
    h+='extern const uint64_t TORSION_PLUS_EVEN_POWER;\n';c+=f'const uint64_t TORSION_PLUS_EVEN_POWER={f};\n'
    for name,ty,val in(('TORSION_ODD_PRIMES','uint64_t',3),('TORSION_ODD_POWERS','uint64_t',b),('TORSION_PLUS_ODD_PRIMES','uint64_t',3),('TORSION_PLUS_ODD_POWERS','size_t',b),('DEGREE_COMMITMENT_POWERS','size_t',b)):
        h+=f'extern const {ty} {name}[1];\n';c+=f'const {ty} {name}[1]={{{val}}};\n'
    h+='extern const ibz_t TORSION_ODD_PRIMEPOWERS[1];\n';c+='const ibz_t TORSION_ODD_PRIMEPOWERS[1]={'+ibz(n)+'};\n'
    for name,v in integer_values.items():h+=f'extern const ibz_t {name};\n';c+=f'const ibz_t {name}={ibz(v)};\n'
    put(inc/'torsion_constants.h',h);put(out/'torsion_constants.c',c)

    O='{'+ibz(2)+','+imat([[2,0,0,1],[0,2,1,0],[0,0,1,0],[0,0,0,1]])+'}'
    qi='{'+ibz(1)+',{'+','.join(ibz(v)for v in(0,1,0,0))+'}}'
    qj='{'+ibz(1)+',{'+','.join(ibz(v)for v in(0,0,1,0))+'}}'
    extremal='{'+O+','+qi+','+qj+',1}'
    put(inc/'quaternion_data.h',banner+'#include <quaternion.h>\n#define NUM_ALTERNATE_EXTREMAL_ORDERS 0\nextern const quat_alg_t QUATALG_PINFTY;\nextern const quat_order_t MAXORD_O0;\nextern const quat_p_extremal_maximal_order_t STANDARD_EXTREMAL_ORDER;\n/* No alternate orders are used or advertised by this implementation. */\nextern const quat_p_extremal_maximal_order_t ALTERNATE_EXTREMAL_ORDERS[1];')
    gram=[[0]*4 for _ in range(4)]
    for k,v in enumerate((1,1,p,p)):gram[k][k]=v
    put(out/'quaternion_data.c',banner+'#include <stdint.h>\n#include <quaternion_data.h>\nconst quat_alg_t QUATALG_PINFTY={'+ibz(p)+','+imat(gram)+'};\nconst quat_order_t MAXORD_O0='+O+';\nconst quat_p_extremal_maximal_order_t STANDARD_EXTREMAL_ORDER='+extremal+';\n/* Storage compatibility only; count is zero. */\nconst quat_p_extremal_maximal_order_t ALTERNATE_EXTREMAL_ORDERS[1]={'+extremal+'};')

    nqr=[];znqr=[];counter=0
    while len(nqr)<20 or len(znqr)<20:
        digest=hashlib.shake_256(f'INKE manuscript field table {level} {counter}'.encode()).digest(16*w)
        z=F2(int.from_bytes(digest[:8*w],'little'),int.from_bytes(digest[8*w:],'little'));counter+=1
        sq=z.sqrt() is not None
        if not sq and len(nqr)<20:nqr.append(z)
        if sq and (z-1).sqrt() is None and len(znqr)<20:znqr.append(z)
    put(inc/'gf_constants.h',banner+'#include <fp2.h>\nextern const fp2_t NQR_TABLE[20];\nextern const fp2_t Z_NQR_TABLE[20];')
    put(out/'gf_constants.c',banner+'#include <stdint.h>\n#include <gf_constants.h>\nconst fp2_t NQR_TABLE[20]={'+','.join(map(fe,nqr))+'};\nconst fp2_t Z_NQR_TABLE[20]={'+','.join(map(fe,znqr))+'};')

    old=out.parent/f'inke_lvl{level}'
    # These parameter-independent trial constants support RepresentInteger.
    # Legacy SQISIGN declarations remain for building shared source files;
    # the manuscript Keygen does not call those legacy algorithms.
    for name in('klpt_constants.c','include/klpt_constants.h'):
        shutil.copyfile(old/name,out/name)
    put(out/'poke_parameters.txt',f'lvl = {level}\np = {hex(p)}\nB = 3\nCfactor = 1')

    def strategy(n,p_cost=20,q_cost=13):
        S={1:[]};C={1:0}
        for i in range(2,n+1):
            j,cost=min(((j,C[i-j]+C[j]+j*p_cost+(i-j)*q_cost)for j in range(1,i)),key=lambda t:t[1])
            S[i]=[j]+S[i-j]+S[j];C[i]=cost
        return S[n]
    def ints(xs):return '{'+','.join(map(str,xs))+'}'
    h=banner+'#ifndef EC_PARAMS_H\n#define EC_PARAMS_H\n#include <stdint.h>\n#include <fp_constants.h>\n#include <inke_field_constants.h>\n'
    defines={'POWER_OF_2':f,'POWER_OF_3':b,'POWER_OF_5':1,'THREEe':40,'FIVEe':27,
             'THREEpF_bitlen':n.bit_length(),'FIVEpF_bitlen':3,'THREE_FIVE_bitlen':(5*n).bit_length(),
             'scaled':1,'gap':83,'P_LEN':1,'M_LEN':0,'sI_max':1,'sJ_max':0,'sK_max':83}
    for name,v in defines.items():h+=f'#define {name} {v}\n'
    arrays={'TWOpF':2**f,'TWOpFm1':2**(f-1),'THREEpE':3**40,'THREEpF':n,'FIVEpE':5**27,'FIVEpF':5,'THREE_FIVE_pF':5*n,'THREEpFdiv2':n//2,'FIVEpFdiv2':2}
    for name,v in arrays.items():h+=f'static digit_t {name}[NWORDS_ORDER]={limbs(v)};\n'
    for name,divisor in(('2f',2**f),('3g',n),('5h',5),('6fg',2**f*n),('35gh',5*n),('235fgh',2**f*n*5)):
        co=(p+1)//divisor;count=max(1,(co.bit_length()+63)//64)
        h+=f'static digit_t p_cofactor_for_{name}[{count}]={limbs(co,count)};\n#define P_COFACTOR_FOR_{name.upper()}_BITLENGTH {co.bit_length()}\n'
    h+='static digit_t p_plus_minus_bitlength[1]={2};\nstatic int sizeI[1]={1};\nstatic int sizeJ[1]={0};\nstatic int sizeK[1]={1};\n'
    h+=f'static int STRATEGY3[{b-1}]={ints(strategy(b,53200,27420))};\n'
    rows=f//2+10
    h+=f'static int STRATEGY4[{rows}][{f//2}]={{\n'+',\n'.join(ints(strategy((f-i)//2))for i in range(rows))+'};\n'
    h+=f'static int strategies[{rows}][{f-1}]={{\n'+',\n'.join(ints(strategy(f-i,16,8))for i in range(rows))+'};\n#endif'
    put(inc/'ec_params.h',h)
    # Keep legacy serialization capacities conservative, with correct field width.
    fpbytes=8*w;fp2bytes=16*w;quatbits=4*p.bit_length();quatbytes=(quatbits+7)//8
    sizes={'FP_ENCODED_BYTES':fpbytes,'FP2_ENCODED_BYTES':fp2bytes,'EC_CURVE_ENCODED_BYTES':fp2bytes,'EC_POINT_ENCODED_BYTES':fp2bytes,'EC_BASIS_ENCODED_BYTES':3*fp2bytes,'CHAIN_LENGTH':8,'QUAT_ALG_ELEM_ENCODED_BITS':quatbits,'QUAT_ALG_ELEM_ENCODED_BYTES':quatbytes,'ID2ISO_LONG_TWO_ISOG_ENCODED_BYTES':8*(2*fp2bytes+2),'ZIP_CHAIN_LEN':12,'ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES':(f+7)//8,'ID2ISO_COMPRESSED_LONG_TWO_ISOG_BYTES':12*((f+7)//8)+1,'SIGNATURE_LEN':12*((f+7)//8)+1+((2**f*n).bit_length()+7)//8+1+(f+8)//8+(n.bit_length()+7)//8,'PUBLICKEY_BYTES':fp2bytes,'SECRETKEY_BYTES':8*fp2bytes+5*quatbytes}
    put(inc/'encoded_sizes.h',banner+'\n'.join(f'#define {k} {v}'for k,v in sizes.items()))
    put(out/'CMakeLists.txt','''# Generated constants for the manuscript uniform-ideal Keygen.
set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF gf_constants.c torsion_constants.c endomorphism_action.c klpt_constants.c quaternion_data.c)
add_library(${LIB_PRECOMP_${SVARIANT_UPPER}} ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${INC_INTBIG} ${INC_QUATERNION} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/ec/ref/include ${PROJECT_SOURCE_DIR}/src/ec/ref/${SVARIANT_LOWER}/include ${PROJECT_SOURCE_DIR}/src/gf/ref/include ${PROJECT_SOURCE_DIR}/src/gf/ref/${SVARIANT_LOWER}/include ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})''')
    put(out/'README.md',f'''These constants implement the manuscript level {level} field, a={a}, b={b},
v2(p+1)={f}, v3(p+1)={data['f3']}. Regenerate with
`python scripts/inke_precompute.py --level {level}` from the implementation root.

Generation is offline only. Runtime arithmetic, uniform sampling, conversion,
Keygen, encryption and decryption use C. The JSON records ordinary affine
coordinates. C field coordinates use Montgomery representation.

The generator checks full torsion orders, basis independence, every evaluated
endomorphism action and the maximal-order relations at the FULL 2-power.
In particular (i+j)/2 is evaluated on a half of the input point in Fp4, never
by arbitrarily halving an output. Alternate extremal orders are not provided;
their advertised count is zero. Legacy SQISIGN tuning/serialization definitions
are retained only to compile shared upstream sources, not as supported INKE
algorithms or wire formats. The uniform Keygen uses the exact O0 and p here.
''')
    print(f'  emitted C constants in {out}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--level',type=int,choices=PARAMS,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--emit-only',action='store_true')
    args=parser.parse_args()
    destination=args.output or Path(__file__).resolve().parents[1]/'src/precomp/ref'/f'inke_uniform_lvl{args.level}'
    start=time.time()
    data=json.loads((destination/'constants.json').read_text()) if args.emit_only else run(args.level,destination)
    emit(data,destination)
    print(f'elapsed {time.time()-start:.3f}s',flush=True)
