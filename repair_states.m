function U = repair_states(U)
%REPAIR_STATES Minimal positivity repair for density and temperature.
%
% This is deliberately conservative and only prevents catastrophic failure.
% A paper-grade code should use a realizability-preserving limiter.

if isvector(U)
    U = U(:);
    vectorInput = true;
else
    vectorInput = false;
end

[m,N] = size(U);
for i = 1:N
    rho = U(1,i);
    if rho < 1e-12 || ~isfinite(rho)
        rho = 1e-12;
        U(:,i) = primitive_to_raw(rho,0,1,m-1);
        continue;
    end

    u = U(2,i)/rho;
    theta = U(3,i)/rho - u^2;

    if theta < 1e-10 || ~isfinite(theta)
        theta = 1e-10;
        U(1,i) = rho;
        U(2,i) = rho*u;
        U(3,i) = rho*(u^2 + theta);
        % Higher moments are left untouched unless nonfinite.
    end

    if any(~isfinite(U(:,i)))
        U(:,i) = primitive_to_raw(rho,u,theta,m-1);
    end
end

if vectorInput
    U = U(:);
end
end
