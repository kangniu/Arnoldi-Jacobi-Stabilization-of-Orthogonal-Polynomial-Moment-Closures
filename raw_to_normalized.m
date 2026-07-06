function W = raw_to_normalized(U,nMax)
%RAW_TO_NORMALIZED Central normalized moments W_0,...,W_nMax.
%
% U can be vector or m-by-N matrix.  Output W is (nMax+1)-by-N.

if isvector(U)
    U = U(:);
end

[m,N] = size(U);
if m < nMax+1
    error('Need raw moments up to U_%d.',nMax);
end

[rho,u,theta] = raw_to_primitive(U);
W = zeros(nMax+1,N);

for k = 0:nMax
    val = zeros(1,N);
    for j = 0:k
        val = val + nchoosek(k,j)*(-u).^(k-j).*U(j+1,:);
    end
    W(k+1,:) = val ./ rho ./ theta.^(k/2);
end

% Enforce exact normalized invariants for numerical cleanliness.
if nMax >= 0, W(1,:) = 1; end
if nMax >= 1, W(2,:) = 0; end
if nMax >= 2, W(3,:) = 1; end

end
