function [rho,u,theta] = raw_to_primitive(U)
%RAW_TO_PRIMITIVE Convert raw moments to rho,u,theta.
%
% U may be a vector or an m-by-N matrix.

rho = U(1,:);
rho = max(rho,1e-12);

u = U(2,:)./rho;
theta = U(3,:)./rho - u.^2;
theta = max(theta,1e-10);

end
