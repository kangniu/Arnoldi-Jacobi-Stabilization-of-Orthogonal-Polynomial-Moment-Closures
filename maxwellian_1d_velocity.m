function M = maxwellian_1d_velocity(v,rho,u,theta)
%MAXWELLIAN_1D_VELOCITY One-dimensional Maxwellian density in velocity.
%
% v     : Nv-by-1 velocity grid
% rho,u,theta : scalars or 1-by-Nx row vectors
%
% M is Nv-by-Nx and satisfies approximately
%   int M dv = rho,
%   int v M dv = rho u,
%   int v^2 M dv = rho (u^2+theta)
% on a sufficiently wide velocity grid.

v = v(:);
rho = rho(:).';
u = u(:).';
theta = max(theta(:).',1e-14);

M = rho ./ sqrt(2*pi*theta) .* exp(-0.5*((v-u).^2 ./ theta));
end
