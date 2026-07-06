function cAsc = b8_charpoly_coeffs_normalized_fd(W3,W4,W5,W6,W7)
%B8_CHARPOLY_COEFFS_NORMALIZED_FD Legacy finite-difference B8 polynomial.
%
% This routine is kept only as a diagnostic comparator.  The production B8
% Jacobi path uses b8_charpoly_coeffs_normalized, which propagates analytic
% sensitivities through the recurrence-based closing moment.

U = [1;0;1;W3;W4;W5;W6;W7];
n = numel(U);
d = zeros(n,1);

epsFD = 2e-6;
for k = 1:n
    h = epsFD*max(1,abs(U(k)));
    Up = U; Um = U;
    Up(k) = Up(k) + h;
    Um(k) = Um(k) - h;
    Up = repair_states(Up);
    Um = repair_states(Um);
    fp = closure_raw_moment(Up,'B8');
    fm = closure_raw_moment(Um,'B8');
    d(k) = (fp-fm)/(Up(k)-Um(k));
end

cAsc = [-d(:); 1];
end
