function [aTail,bTail,info] = recover_jacobi_tail_from_split(Pfull,Pm,Pm1,bm,tailDegree)
%RECOVER_JACOBI_TAIL_FROM_SPLIT Recover a Jacobi tail from a split determinant.
%
% For a Jacobi matrix split after P_m,
%   P_N = A_L P_m - b_m B_{L-1} P_{m-1},
% where L=N-m.  This routine first recovers B_{L-1} from the remainder
% modulo P_m and then applies the Euclidean recurrence to recover the tail
% coefficients of A_L.

Pfull = Pfull(:);
Pm = Pm(:);
Pm1 = Pm1(:);
L = tailDegree;

[~,rFull] = poly_divide_ascending(Pfull,Pm);
rFull = poly_pad_ascending(rFull,numel(Pm)-1);

Rcols = zeros(numel(Pm)-1,L);
for j = 0:L-1
    shifted = [zeros(j,1); Pm1(:)];
    [~,rj] = poly_divide_ascending(shifted,Pm);
    Rcols(:,j+1) = poly_pad_ascending(rj,numel(Pm)-1);
end

% B_{L-1}=lambda^{L-1}+beta_{L-2} lambda^{L-2}+...+beta_0.
rhs = -rFull/bm - Rcols(:,L);
M = Rcols(:,1:L-1);
betaLow = M\rhs;
BTail = [betaLow(:); 1];

term = bm*conv(BTail,Pm1);
Ppad = poly_pad_ascending(Pfull,max(numel(Pfull),numel(term)));
Tpad = poly_pad_ascending(term,numel(Ppad));
[ATail,remA] = poly_divide_ascending(Ppad + Tpad,Pm);
ATail = poly_pad_ascending(ATail,L+1);

[aTail,bTail,chainResidual] = recover_tridiag_from_adjacent_polys(ATail,BTail);

info = struct();
info.ATail = ATail;
info.BTail = BTail;
info.remainderFull = norm(M*betaLow-rhs,2);
info.remainderA = norm(remA,2);
info.chainResidual = chainResidual;
end

function [a,b,residual] = recover_tridiag_from_adjacent_polys(A,B)
% Recover recurrence coefficients from consecutive monic determinants.

A = poly_trim_ascending(A);
B = poly_trim_ascending(B);
L = numel(A)-1;

a = zeros(L,1);
b = zeros(max(L-1,0),1);
residual = 0;

Pcur = A(:);
Pnext = B(:);

for k = 1:L
    if numel(Pnext) == 1
        q = Pcur/Pnext;
        q = poly_pad_ascending(q,2);
        a(k) = -q(1);
        break;
    end

    [q,r] = poly_divide_ascending(Pcur,Pnext);
    q = poly_pad_ascending(q,2);
    a(k) = -q(1);

    r = poly_pad_ascending(r,numel(Pnext)-1);
    coeff = -r(end);
    if coeff <= 1e-12
        coeff = max(abs(coeff),1e-12);
    end
    b(k) = coeff;
    Pafter = r/(-coeff);
    qr = conv(q,Pnext);
    npad = max([numel(Pcur),numel(qr),numel(r)]);
    residual = max(residual,norm(poly_pad_ascending(Pcur,npad) - ...
        poly_pad_ascending(qr,npad) - poly_pad_ascending(r,npad),2));

    Pcur = Pnext;
    Pnext = Pafter;
end
end
