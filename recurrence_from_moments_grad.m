function [a,da,b,db,polys,dpolys,norms,dnorms] = recurrence_from_moments_grad(W,dW,N)
%RECURRENCE_FROM_MOMENTS_GRAD Stieltjes recurrence with sensitivities.
%
% This is the differentiated counterpart of recurrence_from_moments.  It is
% intentionally small and explicit because it is used to remove finite
% differencing from the B8 characteristic-polynomial construction.

W = W(:);
if numel(W) < 2*N
    error('Need moments through W_%d to build %d recurrence coefficients.',2*N-1,N);
end

nvar = size(dW,2);
a = zeros(N,1);
da = zeros(N,nvar);
b = zeros(max(N-1,0),1);
db = zeros(max(N-1,0),nvar);
polys = cell(N+1,1);
dpolys = cell(N+1,1);
norms = zeros(N+1,1);
dnorms = zeros(N+1,nvar);

polys{1} = 1;
dpolys{1} = zeros(1,nvar);
[norms(1),dnorms(1,:)] = moment_poly_inner_grad(polys{1},dpolys{1},polys{1},dpolys{1},W,dW);

if N == 0
    return;
end

pm = 0;
dpm = zeros(1,nvar);
pk = polys{1};
dpk = dpolys{1};

for k = 0:N-1
    [normk,dnormk] = moment_poly_inner_grad(pk,dpk,pk,dpk,W,dW);
    norms(k+1) = normk;
    dnorms(k+1,:) = dnormk;

    xpk = [0; pk(:)];
    dxpk = [zeros(1,nvar); dpk];
    [num, dnum] = moment_poly_inner_grad(xpk,dxpk,pk,dpk,W,dW);
    a(k+1) = num/normk;
    da(k+1,:) = (dnum*normk - num*dnormk)/(normk^2);

    pnext = [-a(k+1)*pk(:); 0] + [0; pk(:)];
    dpnext = [-da(k+1,:).*pk(:) - a(k+1)*dpk; zeros(1,nvar)] + ...
        [zeros(1,nvar); dpk];

    if k > 0
        [old,dold] = pad_poly_grad(pm,dpm,numel(pnext));
        pnext = pnext - b(k)*old;
        dpnext = dpnext - db(k,:).*old - b(k)*dold;
    end

    polys{k+2} = pnext(:);
    dpolys{k+2} = dpnext;

    if k < N-1
        [normNext,dnormNext] = moment_poly_inner_grad(pnext,dpnext,pnext,dpnext,W,dW);
        norms(k+2) = normNext;
        dnorms(k+2,:) = dnormNext;
        b(k+1) = normNext/normk;
        db(k+1,:) = (dnormNext*normk - normNext*dnormk)/(normk^2);
    end

    pm = pk;
    dpm = dpk;
    pk = pnext;
    dpk = dpnext;
end

[norms(N+1),dnorms(N+1,:)] = moment_poly_inner_grad(polys{N+1},dpolys{N+1},polys{N+1},dpolys{N+1},W,dW);
end

function [pout,dpout] = pad_poly_grad(p,dp,n)
p = p(:);
nvar = size(dp,2);
pout = zeros(n,1);
dpout = zeros(n,nvar);
m = min(numel(p),n);
pout(1:m) = p(1:m);
dpout(1:m,:) = dp(1:m,:);
end
