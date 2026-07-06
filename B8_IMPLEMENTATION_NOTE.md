# B8 implementation note

This package implements a complete B4/B6/B8 V+A/Arnoldi--Jacobi workflow.

The JCP paper gives the B8 closure-coefficient pattern but does not print the
resulting closing flux or tail recurrence coefficients because the symbolic
expressions are very large.  The implementation here therefore does not
hand-transcribe a B8 formula.  Instead:

1. `b8_closing_moment_normalized.m` constructs the B8 closing moment from
   the recurrence definition of the hierarchy coefficient `b4`.
2. `b8_closing_moment_normalized.m` also propagates forward sensitivities of
   the recurrence construction, giving `dW8/d[W3,...,W7]`.
3. `b8_charpoly_coeffs_normalized.m` uses these sensitivities and the
   normalized-to-raw moment chain rule to form the B8 characteristic
   polynomial coefficients directly.
4. `b8_jacobi_coefficients.m` recovers the Jacobi tail from the split
   polynomial identity `P8 = A4*P4 - b4*B3*P3`.

The B8 implementation is validated by the same tests used for B6:

1. equilibrium recovery,
2. positive Jacobi coefficients,
3. static Jacobian-vs-Jacobi wave-speed agreement,
4. small-grid PDE agreement,
5. grid/time sweep.

The legacy finite-difference polynomial generator is kept as
`b8_charpoly_coeffs_normalized_fd.m` only for diagnostics.  The production
V+A/Jacobi B8 path does not use centered finite differences to construct
`P8`.
