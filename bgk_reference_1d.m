function ref = bgk_reference_1d(cfg,varargin)
%BGK_REFERENCE_1D Discrete-velocity BGK reference solver for the Riemann case.
%
% This is a compact reference implementation for article diagnostics.  It is
% not intended to replace a high-resolution kinetic production code, but it
% provides an independent kinetic reference for the moment-closure profiles.
%
% Name-value options:
%   'Nx'       : spatial cells, default cfg.Nx
%   'Nv'       : velocity cells, default 501
%   'vmin'     : minimum velocity, default -5000
%   'vmax'     : maximum velocity, default  5000
%   'tEnd'     : final time, default cfg.tEnd
%   'tau'      : BGK relaxation time, default cfg.tau
%   'CFL'      : CFL number, default min(cfg.CFL,0.8)
%   'printEvery': progress frequency, default 200
%
% Output fields include x, v, f, U, rho, u, theta, W, status, nSteps, cpu.

opts = struct();
opts.Nx = cfg.Nx;
opts.Nv = 501;
opts.vmin = -5000;
opts.vmax = 5000;
opts.tEnd = cfg.tEnd;
opts.tau = cfg.tau;
opts.CFL = min(cfg.CFL,0.8);
opts.printEvery = 200;

if mod(numel(varargin),2) ~= 0
    error('Options must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    val = varargin{k+1};
    if ~isfield(opts,name)
        error('Unknown option: %s',name);
    end
    opts.(name) = val;
end

xmin = cfg.xmin; xmax = cfg.xmax;
Nx = opts.Nx; Nv = opts.Nv;
x = linspace(xmin,xmax,Nx);
dx = (xmax-xmin)/(Nx-1);

v = linspace(opts.vmin,opts.vmax,Nv).';
wv = velocity_quadrature_weights(v);

thetaL = cfg.left.p/cfg.left.rho;
thetaR = cfg.right.p/cfg.right.rho;

fL = maxwellian_1d_velocity(v,cfg.left.rho,cfg.left.u,thetaL);
fR = maxwellian_1d_velocity(v,cfg.right.rho,cfg.right.u,thetaR);

f = zeros(Nv,Nx);
for i = 1:Nx
    if x(i) <= 0
        f(:,i) = fL;
    else
        f(:,i) = fR;
    end
end

vmaxAbs = max(abs(v));
dtAdv = opts.CFL*dx/max(vmaxAbs,eps);
t = 0;
nSteps = 0;
status = 'finished';

pos = find(v > 0);
neg = find(v < 0);
vp = v(pos);
vn = v(neg);

fprintf('\nBGK discrete-velocity reference\n');
fprintf('  Nx=%d, Nv=%d, tEnd=%.6e, tau=%g, v=[%.1f, %.1f]\n', ...
    Nx,Nv,opts.tEnd,opts.tau,opts.vmin,opts.vmax);

tic;
while t < opts.tEnd - 1e-15
    dt = min(dtAdv,opts.tEnd-t);
    fOld = f;
    fNew = fOld;

    % Positive velocities: upwind from left.
    if ~isempty(pos)
        fNew(pos,1) = fOld(pos,1) - (dt/dx) .* (vp .* fOld(pos,1) - vp .* fL(pos));
        if Nx > 1
            fNew(pos,2:Nx) = fOld(pos,2:Nx) - (dt/dx) .* ...
                (vp .* fOld(pos,2:Nx) - vp .* fOld(pos,1:Nx-1));
        end
    end

    % Negative velocities: upwind from right.
    if ~isempty(neg)
        if Nx > 1
            fNew(neg,1:Nx-1) = fOld(neg,1:Nx-1) - (dt/dx) .* ...
                (vn .* fOld(neg,2:Nx) - vn .* fOld(neg,1:Nx-1));
        end
        fNew(neg,Nx) = fOld(neg,Nx) - (dt/dx) .* (vn .* fR(neg) - vn .* fOld(neg,Nx));
    end

    fNew = max(fNew,0);

    % Point-implicit BGK relaxation.
    if isfinite(opts.tau) && opts.tau > 0
        U2 = moments_from_distribution(v,wv,fNew,2);
        rho = max(U2(1,:),1e-14);
        u = U2(2,:)./rho;
        theta = max(U2(3,:)./rho - u.^2,1e-14);
        M = maxwellian_1d_velocity(v,rho,u,theta);
        alpha = dt/opts.tau;
        fNew = (fNew + alpha*M)/(1+alpha);
    end

    f = fNew;
    t = t + dt;
    nSteps = nSteps + 1;

    if opts.printEvery > 0 && mod(nSteps,opts.printEvery)==0
        U2 = moments_from_distribution(v,wv,f,2);
        rhoTmp = U2(1,:);
        thetaTmp = U2(3,:)./rhoTmp - (U2(2,:)./rhoTmp).^2;
        fprintf('      BGK step=%5d, t=%.4e, dt=%.2e, minrho=%.3e, mintheta=%.3e\n', ...
            nSteps,t,dt,min(rhoTmp),min(thetaTmp));
    end

    if nSteps > 200000
        status = 'aborted_max_steps';
        break;
    end
end
cpu = toc;

U = moments_from_distribution(v,wv,f,5);
[rho,u,theta] = raw_to_primitive(U(1:3,:));
W = raw_to_normalized(U,5);

ref = struct();
ref.x = x(:);
ref.v = v;
ref.wv = wv;
ref.f = f;
ref.U = U;
ref.rho = rho;
ref.u = u;
ref.theta = theta;
ref.W = W;
ref.t = t;
ref.nSteps = nSteps;
ref.status = status;
ref.cpu = cpu;
ref.options = opts;

fprintf('  BGK status=%s, steps=%d, final time %.6e, cpu %.2fs\n', ...
    status,nSteps,t,cpu);
fprintf('  BGK min rho %.4e, min theta %.4e\n',min(rho),min(theta));

end
