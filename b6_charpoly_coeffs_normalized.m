function cAsc = b6_charpoly_coeffs_normalized(W3,W4,W5)
%B6_CHARPOLY_COEFFS_NORMALIZED Characteristic polynomial coefficients for B6.
%
% For normalized state rho=1, u=0, theta=1, the B6 flux Jacobian has
% companion form.  If d_k = dU6/dU_k at the normalized state, then
%
%   p(lambda) = lambda^6 - d5 lambda^5 - ... - d1 lambda - d0.
%
% This routine returns coefficients in ascending order:
%   cAsc = [-d0,-d1,-d2,-d3,-d4,-d5,1].
%
% The derivatives are generated from the B6 closing flux obtained by solving
% the recurrence condition
%
%   b3 = (b1+b2) + 1/2[(a2-a0)^2+(a2-a1)^2].
%
% This corrected form gives W6=15 at equilibrium.

w3 = W3; w4 = W4; w5 = W5;
g = w3^2 - w4 + 1;

if abs(g) < 1e-12
    g = sign_nonzero(g)*1e-12;
end

d0 = -0.5*(5*w3^6 - 20*w3^4*w4 + 2*w3^4 + ...
    8*w3^3*w5 + 21*w3^2*w4^2 - 6*w3^2*w4 + w3^2 - ...
    16*w3*w4*w5 - 2*w4^3 + 4*w4^2 - 2*w4 + 4*w5^2)/g^2;

d1 = (4*w3^7 - 14*w3^5*w4 + 13*w3^5 - w3^4*w5 + ...
    20*w3^3*w4^2 - 34*w3^3*w4 + 6*w3^3 + ...
    2*w3^2*w4*w5 + 10*w3^2*w5 - 14*w3*w4^3 + ...
    21*w3*w4^2 - 8*w3*w4 - 4*w3*w5^2 + w3 + ...
    7*w4^2*w5 - 10*w4*w5 + 3*w5)/g^2;

d2 = 0.5*(9*w3^8 - 32*w3^6*w4 + 18*w3^6 + ...
    6*w3^5*w5 + 37*w3^4*w4^2 - 56*w3^4*w4 - 3*w3^4 - ...
    12*w3^3*w4*w5 + 28*w3^3*w5 - 14*w3^2*w4^3 + ...
    40*w3^2*w4^2 + 6*w3^2*w4 - 4*w3^2*w5^2 + ...
    14*w3*w4^2*w5 - 44*w3*w4*w5 - 2*w3*w5 - ...
    4*w4^4 + 6*w4^3 - 2*w4 + 8*w5^2)/g^2;

d3 = -(6*w3^7 - 19*w3^5*w4 + 7*w3^5 + 3*w3^4*w5 + ...
    20*w3^3*w4^2 - 16*w3^3*w4 - 4*w3^3 - ...
    2*w3^2*w4*w5 + 10*w3^2*w5 - 11*w3*w4^3 + ...
    5*w3*w4^2 + 7*w3*w4 - 4*w3*w5^2 - w3 + ...
    7*w4^2*w5 - 6*w4*w5 - w5)/g^2;

d4 = 0.5*(7*w3^6 - 14*w3^4*w4 + 4*w3^4 + ...
    8*w3^3*w5 - 5*w3^2*w4^2 - 2*w3^2*w4 - 9*w3^2 + ...
    16*w3*w5 + 8*w4^3 - 18*w4^2 + 12*w4 - ...
    4*w5^2 - 2)/g^2;

d5 = -(3*w3^3 - 7*w3*w4 - w3 + 4*w5)/g;

cAsc = [-d0; -d1; -d2; -d3; -d4; -d5; 1];

end

function s = sign_nonzero(x)
if x >= 0
    s = 1;
else
    s = -1;
end
end
