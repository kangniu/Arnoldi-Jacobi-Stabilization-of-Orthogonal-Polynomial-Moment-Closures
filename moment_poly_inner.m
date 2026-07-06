function val = moment_poly_inner(p,q,W)
%MOMENT_POLY_INNER Inner product of two polynomials from moments.
%
% p, q, and W use ascending powers. W(k+1) is the moment W_k.

p = p(:);
q = q(:);
W = W(:);

val = 0;
for i = 1:numel(p)
    for j = 1:numel(q)
        idx = (i-1) + (j-1) + 1;
        if idx > numel(W)
            error('Need moment W_%d for polynomial inner product.',idx-1);
        end
        val = val + p(i)*q(j)*W(idx);
    end
end
end
