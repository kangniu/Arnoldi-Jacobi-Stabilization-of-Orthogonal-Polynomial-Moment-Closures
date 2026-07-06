function lambda = wave_speeds_B4_jacobi(U)
%WAVE_SPEEDS_B4_JACOBI Jacobi/Arnoldi wave speeds for B4 closure.
%
% In normalized variables:
%   a0=0, b1=1,
%   a1=W3, b2=W4-W3^2-1 = W3^2+2 for B4,
%   a2=W3(3W3^2+4)/(W3^2+2),
%   a3=2W3/(W3^2+2),
%   b3=(9W3^4+20W3^2+12)/(W3^2+2)^2.
%
% The physical speeds are u + sqrt(theta)*eig(J).

[rho,u,theta] = raw_to_primitive(U);
W = raw_to_normalized(U,3);
W3 = W(4);

a0 = 0;
b1 = 1;
a1 = W3;
b2 = W3^2 + 2;
a2 = W3*(3*W3^2 + 4)/(W3^2 + 2);
a3 = 2*W3/(W3^2 + 2);
b3 = (9*W3^4 + 20*W3^2 + 12)/(W3^2 + 2)^2;

J = diag([a0,a1,a2,a3]) + ...
    diag(sqrt([b1,b2,b3]),1) + diag(sqrt([b1,b2,b3]),-1);

lambdaHat = eig((J+J')/2);
lambda = u + sqrt(theta)*lambdaHat;
lambda = sort(real(lambda));
end
