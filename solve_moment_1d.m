function sol = solve_moment_1d(U,x,tEnd,tau,closure,CFL,waveMethod,stop)
%SOLVE_MOMENT_1D First-order finite-volume HLL solver with implicit BGK source.
%
% Stopping criteria fields in STOP:
%   maxSteps, minDt, maxWaveSpeed, abortOnSmallDt, abortOnLargeWaveSpeed,
%   printEvery.

if nargin < 7 || isempty(waveMethod)
    waveMethod = 'jacobian';
end
if nargin < 8 || isempty(stop)
    stop = struct();
end

if ~isfield(stop,'maxSteps'), stop.maxSteps = 20000; end
if ~isfield(stop,'minDt'), stop.minDt = 1e-12; end
if ~isfield(stop,'maxWaveSpeed'), stop.maxWaveSpeed = Inf; end
if ~isfield(stop,'abortOnSmallDt'), stop.abortOnSmallDt = true; end
if ~isfield(stop,'abortOnLargeWaveSpeed'), stop.abortOnLargeWaveSpeed = true; end
if ~isfield(stop,'printEvery'), stop.printEvery = 50; end

dx = x(2)-x(1);
N = size(U,2);
t = 0;
nSteps = 0;
status = 'finished';

while t < tEnd
    if nSteps >= stop.maxSteps
        status = 'aborted_max_steps';
        warning('Stopping: reached maxSteps=%d at t=%.6e.',stop.maxSteps,t);
        break;
    end

    % Compute stable time step.
    smax = 0;
    badLambda = false;
    for i = 1:N
        lam = wave_speeds(U(:,i),closure,waveMethod);
        if any(~isfinite(lam))
            badLambda = true;
            break;
        end
        smax = max(smax,max(abs(lam)));
    end

    if badLambda || ~isfinite(smax)
        status = 'aborted_bad_wave_speed';
        warning('Stopping: nonfinite wave speed detected at step=%d, t=%.6e.',nSteps,t);
        break;
    end

    if smax < 1e-12
        smax = 1e-12;
    end

    if stop.abortOnLargeWaveSpeed && smax > stop.maxWaveSpeed
        status = 'aborted_large_wave_speed';
        warning('Stopping: max wave speed %.4e exceeds threshold %.4e at step=%d, t=%.6e.', ...
            smax,stop.maxWaveSpeed,nSteps,t);
        break;
    end

    dt = CFL*dx/smax;
    if t + dt > tEnd
        dt = tEnd - t;
    end

    if stop.abortOnSmallDt && dt < stop.minDt
        status = 'aborted_small_dt';
        warning('Stopping: dt %.4e below minDt %.4e at step=%d, t=%.6e, maxs=%.4e.', ...
            dt,stop.minDt,nSteps,t,smax);
        break;
    end

    % Ghost cells: outflow boundary.
    Ug = [U(:,1), U, U(:,end)];
    Fh = zeros(size(U,1),N+1);

    for i = 1:N+1
        UL = Ug(:,i);
        UR = Ug(:,i+1);
        Fh(:,i) = hll_flux(UL,UR,closure,waveMethod);
    end

    Ustar = U;
    for i = 1:N
        Ustar(:,i) = U(:,i) - dt/dx*(Fh(:,i+1)-Fh(:,i));
    end

    Ustar = repair_states(Ustar);

    % Point-implicit BGK relaxation for non-conserved moments.
    U = bgk_relax_implicit(Ustar,dt,tau);
    U = repair_states(U);

    t = t + dt;
    nSteps = nSteps + 1;

    if stop.printEvery > 0 && mod(nSteps,stop.printEvery)==0
        [rho,~,theta] = raw_to_primitive(U);
        fprintf('      step=%5d, t=%.4e, dt=%.2e, maxs=%.2e, minrho=%.2e, mintheta=%.2e\n', ...
            nSteps,t,dt,smax,min(rho),min(theta));
    end
end

sol.x = x;
sol.U = U;
sol.t = t;
sol.nSteps = nSteps;
sol.closure = closure;
sol.waveMethod = waveMethod;
sol.status = status;

end
