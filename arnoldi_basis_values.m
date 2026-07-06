function [Q,H] = arnoldi_basis_values(x,n,reorth)
%ARNOLDI_BASIS_VALUES Arnoldi basis for the columns of a Vandermonde matrix.
%
% Q(:,1),...,Q(:,n+1) are orthogonal with respect to the discrete inner
% product <a,b> = a'*b/m.  H stores the Arnoldi recurrence coefficients.
%
% This is a corrected and isolated version of polyfitA's basis-generation
% part.  The reorthogonalization loop is over all previous basis vectors.

if nargin < 3
    reorth = true;
end

x = x(:);
m = numel(x);
Q = zeros(m,n+1);
Q(:,1) = 1;
H = zeros(n+1,n);

for k = 1:n
    q = x .* Q(:,k);

    for j = 1:k
        H(j,k) = Q(:,j)'*q/m;
        q = q - H(j,k)*Q(:,j);
    end

    if reorth
        for j = 1:k
            hcorr = Q(:,j)'*q/m;
            q = q - hcorr*Q(:,j);
            H(j,k) = H(j,k) + hcorr;
        end
    end

    H(k+1,k) = norm(q)/sqrt(m);
    if H(k+1,k) < 1e-14
        error('Arnoldi breakdown at k=%d. Try lower degree or different nodes.',k);
    end
    Q(:,k+1) = q/H(k+1,k);
end

end
