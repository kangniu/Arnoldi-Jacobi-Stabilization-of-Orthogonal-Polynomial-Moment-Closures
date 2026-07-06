function U = moments_from_distribution(v,wv,f,nMax)
%MOMENTS_FROM_DISTRIBUTION Raw velocity moments from a discrete distribution.
%
% v  : Nv-by-1 velocity grid
% wv : Nv-by-1 quadrature weights
% f  : Nv-by-Nx distribution
% U  : (nMax+1)-by-Nx raw moment array

v = v(:);
wv = wv(:);
Nx = size(f,2);
U = zeros(nMax+1,Nx);

vw = wv;
pow = ones(size(v));
for k = 0:nMax
    U(k+1,:) = (vw .* pow).' * f;
    pow = pow .* v;
end
end
