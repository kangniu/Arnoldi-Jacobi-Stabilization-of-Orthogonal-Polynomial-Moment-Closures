function [a,b,polys,norms] = recurrence_from_moments(W,N)
%RECURRENCE_FROM_MOMENTS Stieltjes recurrence from moments W_0,... .
%
% Returns a=[a0,...,a_{N-1}], b=[b1,...,b_{N-1}], the monic polynomials
% P0,...,PN, and their squared norms.  Coefficients are ascending.

W = W(:);
if numel(W) < 2*N
    error('Need moments through W_%d to build %d recurrence coefficients.',2*N-1,N);
end

a = zeros(N,1);
b = zeros(max(N-1,0),1);
polys = cell(N+1,1);
norms = zeros(N+1,1);

polys{1} = 1;
norms(1) = moment_poly_inner(polys{1},polys{1},W);

if N == 0
    return;
end

pm = 0;
pk = polys{1};

for k = 0:N-1
    normk = moment_poly_inner(pk,pk,W);
    norms(k+1) = normk;
    xpk = [0; pk(:)];
    a(k+1) = moment_poly_inner(xpk,pk,W)/normk;

    pnext = [-a(k+1)*pk(:); 0] + [0; pk(:)];
    if k > 0
        old = poly_pad_ascending(pm,numel(pnext));
        pnext = pnext - b(k)*old;
    end

    polys{k+2} = pnext(:);
    if k < N-1
        normNext = moment_poly_inner(pnext,pnext,W);
        norms(k+2) = normNext;
        b(k+1) = normNext/normk;
    end

    pm = pk;
    pk = pnext;
end

norms(N+1) = moment_poly_inner(polys{N+1},polys{N+1},W);
end
