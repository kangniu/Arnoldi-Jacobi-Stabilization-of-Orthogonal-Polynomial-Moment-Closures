function Uk = central_normalized_to_raw(rho,u,theta,W,k)
%CENTRAL_NORMALIZED_TO_RAW Convert W_0,...,W_k to raw U_k.

if numel(W) < k+1
    error('Need W up to W_%d.',k);
end

val = 0;
for j = 0:k
    val = val + nchoosek(k,j)*u^(k-j)*theta^(j/2)*W(j+1);
end
Uk = rho*val;
end
