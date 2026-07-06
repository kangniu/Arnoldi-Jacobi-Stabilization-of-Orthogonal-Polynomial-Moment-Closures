function p = poly_pad_ascending(p,n)
%POLY_PAD_ASCENDING Pad ascending polynomial to length n.
p = p(:);
if numel(p) < n
    p = [p; zeros(n-numel(p),1)];
elseif numel(p) > n
    p = p(1:n);
end
end
