function [q,r] = poly_divide_ascending(p,d)
%POLY_DIVIDE_ASCENDING Polynomial division for ascending coefficient vectors.
%
% p(lambda) = q(lambda)*d(lambda) + r(lambda), deg(r)<deg(d).

p = poly_trim_ascending(p);
d = poly_trim_ascending(d);

if abs(d(end)) < eps
    error('Divisor leading coefficient is zero.');
end

np = numel(p)-1;
nd = numel(d)-1;

if np < nd
    q = 0;
    r = p;
    return;
end

r = p;
q = zeros(np-nd+1,1);

while numel(r)-1 >= nd
    degR = numel(r)-1;
    coeff = r(end)/d(end);
    shift = degR-nd;
    q(shift+1) = coeff;

    sub = [zeros(shift,1); coeff*d(:)];
    r = r - poly_pad_ascending(sub,numel(r));
    r = poly_trim_ascending(r,1e-12);
end

r = r(:);
q = q(:);
end
