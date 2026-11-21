These constants implement the manuscript level 3 field, a=384, b=243,
v2(p+1)=386, v3(p+1)=244. Regenerate with
`python scripts/inke_precompute.py --level 3` from the implementation root.

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

