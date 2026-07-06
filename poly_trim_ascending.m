function p = poly_trim_ascending(p,tol)
%POLY_TRIM_ASCENDING Trim trailing tiny coefficients in ascending polynomial.
if nargin < 2
    tol = 1e-14;
end
p = p(:);
while numel(p) > 1 && abs(p(end)) <= tol
    p(end) = [];
end
end
