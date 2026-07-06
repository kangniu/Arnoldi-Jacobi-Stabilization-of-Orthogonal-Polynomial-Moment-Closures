function [lambda,info] = wave_speeds_B8_jacobi(U)
%WAVE_SPEEDS_B8_JACOBI Arnoldi/Jacobi wave speeds for the B8 closure.

[~,u,theta] = raw_to_primitive(U);
W = raw_to_normalized(U,7);
W3 = W(4); W4 = W(5); W5 = W(6); W6 = W(7); W7 = W(8);

[a,b,info] = b8_jacobi_coefficients(W3,W4,W5,W6,W7);

J = diag(a) + diag(sqrt(b),1) + diag(sqrt(b),-1);
lambdaHat = eig((J+J')/2);
lambda = sort(real(u + sqrt(theta)*lambdaHat));
end
