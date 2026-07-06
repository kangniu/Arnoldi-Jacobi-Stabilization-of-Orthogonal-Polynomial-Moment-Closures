function U = initialize_riemann_state(x,cfg,nMom)
%INITIALIZE_RIEMANN_STATE Maxwellian left/right raw moments.
N = numel(x);
U = zeros(nMom,N);

thetaL = cfg.left.p/cfg.left.rho;
thetaR = cfg.right.p/cfg.right.rho;

UL = primitive_to_raw(cfg.left.rho,cfg.left.u,thetaL,nMom-1);
UR = primitive_to_raw(cfg.right.rho,cfg.right.u,thetaR,nMom-1);

for i = 1:N
    if x(i) <= 0
        U(:,i) = UL;
    else
        U(:,i) = UR;
    end
end
end
