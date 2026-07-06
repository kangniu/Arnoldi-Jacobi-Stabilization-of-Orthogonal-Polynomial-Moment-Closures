function p = recurrence_poly_ascending(a,b)
%RECURRENCE_POLY_ASCENDING Build monic orthogonal polynomial from recurrence.
%
% a = [a0,...,a_{N-1}], b = [b1,...,b_{N-1}].
% Returns p_N(lambda) in ascending powers:
%   p(1) + p(2) lambda + ... + p(N+1) lambda^N.

N = numel(a);
p0 = 1;              % P0
if N == 0
    p = p0(:);
    return;
end

p1 = [-a(1); 1];     % P1 = lambda - a0
if N == 1
    p = p1(:);
    return;
end

pm = p0(:);
pk = p1(:);

for k = 2:N
    % P_k = (lambda - a_{k-1}) P_{k-1} - b_{k-1} P_{k-2}.
    term = [-a(k)*pk; 0] + [0; pk];
    old = [pm; zeros(numel(term)-numel(pm),1)];
    pnew = term - b(k-1)*old;
    pm = pk;
    pk = pnew;
end

p = pk(:);
end
