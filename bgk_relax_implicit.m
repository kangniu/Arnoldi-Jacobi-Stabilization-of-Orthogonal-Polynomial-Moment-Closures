function Unew = bgk_relax_implicit(U,dt,tau)
%BGK_RELAX_IMPLICIT Point-implicit relaxation toward Maxwellian moments.
%
% Conserved moments U0,U1,U2 are unchanged.  Higher moments are relaxed as
%     U_k^{n+1} = (U_k^* + dt/tau U_k^M)/(1+dt/tau), k>=3.

if isinf(tau) || tau <= 0
    Unew = U;
    return;
end

[m,N] = size(U);
Unew = U;
alpha = dt/tau;

for i = 1:N
    [rho,u,theta] = raw_to_primitive(U(:,i));
    UM = primitive_to_raw(rho,u,theta,m-1);
    for k = 4:m   % MATLAB row k corresponds to moment U_{k-1}; start at U3.
        Unew(k,i) = (U(k,i) + alpha*UM(k))/(1+alpha);
    end
end
end
