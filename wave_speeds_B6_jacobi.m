function [lambda,info] = wave_speeds_B6_jacobi(U)
%WAVE_SPEEDS_B6_JACOBI Prototype Arnoldi/Jacobi wave speeds for B6 closure.
%
% The normalized Jacobi matrix is recovered by matching the B6 characteristic
% polynomial with a three-term recurrence.  Physical speeds are then
%   lambda = u + sqrt(theta) eig(J).
%
% This function is intended for static tests and moderate-size diagnostics.
% It is more expensive than the B4 closed-form Jacobi formula because B6
% coefficient matching is solved numerically.

[rho,u,theta] = raw_to_primitive(U);
W = raw_to_normalized(U,5);
W3 = W(4); W4 = W(5); W5 = W(6);

[a,b,info] = b6_jacobi_coefficients(W3,W4,W5);

J = diag(a) + diag(sqrt(b),1) + diag(sqrt(b),-1);
lambdaHat = eig((J+J')/2);
lambda = sort(real(u + sqrt(theta)*lambdaHat));

end
