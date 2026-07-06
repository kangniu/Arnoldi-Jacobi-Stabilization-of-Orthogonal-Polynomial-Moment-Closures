function J = flux_jacobian_fd(U,closure)
%FLUX_JACOBIAN_FD Finite-difference approximation of dF/dU.
%
% This represents the "original-paper style" computation in this research
% code: build the flux Jacobian and compute its eigenvalues.  The paper's
% flux Jacobian has companion form; here the derivatives of the closing flux
% are approximated numerically.

m = numel(U);
J = zeros(m,m);

epsFD = 1e-6;
for k = 1:m
    h = epsFD*max(1,abs(U(k)));
    Up = U; Um = U;
    Up(k) = Up(k) + h;
    Um(k) = Um(k) - h;

    Up = repair_states(Up);
    Um = repair_states(Um);

    Fp = flux_moment(Up,closure);
    Fm = flux_moment(Um,closure);
    denom = Up(k)-Um(k);
    if abs(denom) < eps
        denom = 2*h;
    end
    J(:,k) = (Fp-Fm)/denom;
end

end
