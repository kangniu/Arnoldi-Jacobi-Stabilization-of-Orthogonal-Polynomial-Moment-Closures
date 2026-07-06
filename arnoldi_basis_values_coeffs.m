function [Q,H,C] = arnoldi_basis_values_coeffs(x,n,reorth)
%ARNOLDI_BASIS_VALUES_COEFFS Arnoldi basis values and monomial coefficients.
%
% Q(:,k+1) = q_k(x_i), where q_k is represented by
%   q_k(t) = C(1,k+1) + C(2,k+1)t + ... + C(n+1,k+1)t^n.
%
% The discrete inner product is <a,b> = a'*b/m, so Q'*Q/m ~= I.

if nargin < 3
    reorth = true;
end

x = x(:);
m = numel(x);
Q = zeros(m,n+1);
Q(:,1) = 1;

H = zeros(n+1,n);
C = zeros(n+1,n+1);
C(1,1) = 1;

for k = 1:n
    q = x .* Q(:,k);

    % coefficient of t*q_{k-1}(t)
    c = zeros(n+1,1);
    c(2:end) = C(1:end-1,k);

    for j = 1:k
        H(j,k) = Q(:,j)'*q/m;
        q = q - H(j,k)*Q(:,j);
        c = c - H(j,k)*C(:,j);
    end

    if reorth
        for j = 1:k
            hcorr = Q(:,j)'*q/m;
            q = q - hcorr*Q(:,j);
            c = c - hcorr*C(:,j);
            H(j,k) = H(j,k) + hcorr;
        end
    end

    H(k+1,k) = norm(q)/sqrt(m);
    if H(k+1,k) < 1e-14
        error('Arnoldi breakdown at k=%d. Try lower degree or different nodes.',k);
    end

    Q(:,k+1) = q/H(k+1,k);
    C(:,k+1) = c/H(k+1,k);
end
end
