function J = legendre_jacobi_matrix(N,center,scale)
%LEGENDRE_JACOBI_MATRIX Orthonormal Legendre Jacobi matrix.
%
% The unshifted matrix has eigenvalues in (-1,1).  The shifted/scaled version
% is center*I + scale*J0, so eigenvalues lie in (center-scale, center+scale).
%
% This provides a controlled model for high-order characteristic-polynomial
% root computations: eig(J) is the stable Arnoldi/Jacobi route, while
% roots(poly(J)) mimics monomial companion/root computation.

if nargin < 2, center = 0; end
if nargin < 3, scale = 1; end

k = (1:N-1)';
beta = k ./ sqrt(4*k.^2 - 1);
J0 = diag(beta,1) + diag(beta,-1);
J = center*eye(N) + scale*J0;
end
