function [a,b,info] = b6_jacobi_coefficients(W3,W4,W5)
%B6_JACOBI_COEFFICIENTS Semi-explicit B6 Jacobi recurrence.
%
% This routine replaces the former fminsearch coefficient matching by a
% recurrence-based semi-explicit construction.
%
% Known recurrence coefficients:
%   a0 = 0, b1 = 1,
%   a1 = W3,
%   b2 = W4 - W3^2 - 1,
%   a2 = (W3^3 - 2 W3 W4 + W5)/b2,
%   b3 = (b1+b2) + 1/2[(a2-a0)^2+(a2-a1)^2].
%
% Let P2 and P3 be the known leading principal characteristic polynomials.
% For the unknown tail coefficients a3,a4,a5,b4,b5, the sixth polynomial
% can be written in block-recursive form
%
%   P6 = A3 P3 - b3 B2 P2,
%
% where
%   B2 = (lambda-a5)(lambda-a4) - b5,
%   A3 = (lambda-a3)B2 - b4(lambda-a5).
%
% Therefore B2 is recovered from the remainder of P6 modulo P3 by solving a
% small 3-by-2 linear problem.  Then A3 is obtained by polynomial division,
% and a3,b4,a5,a4,b5 follow from coefficient comparisons.
%
% This is the intended "semi-explicit recurrence" implementation and avoids
% nonlinear optimization inside B6 wave-speed evaluations.

D = W4 - W3^2 - 1;
if D <= 1e-12
    D = 1e-12;
end

a0 = 0;
b1 = 1;
a1 = W3;
b2 = D;
a2 = (W3^3 - 2*W3*W4 + W5)/D;
b3 = (b1 + b2) + 0.5*((a2-a0)^2 + (a2-a1)^2);
b3 = max(b3,1e-12);

p2 = recurrence_poly_ascending([a0;a1],[b1]);
p3 = recurrence_poly_ascending([a0;a1;a2],[b1;b2]);
p6 = b6_charpoly_coeffs_normalized(W3,W4,W5);

% Remainder of P6 modulo P3.
[~,r6] = poly_divide_ascending(p6,p3);
r6 = poly_pad_ascending(r6,3);

% Remainders of lambda^k P2 modulo P3, k=0,1,2.
[~,R0] = poly_divide_ascending(p2,p3);
[~,R1] = poly_divide_ascending([0;p2],p3);
[~,R2] = poly_divide_ascending([0;0;p2],p3);

R0 = poly_pad_ascending(R0,3);
R1 = poly_pad_ascending(R1,3);
R2 = poly_pad_ascending(R2,3);

% B2 = lambda^2 + beta1 lambda + beta0.
% r6 = -b3*(R2 + beta1*R1 + beta0*R0).
rhs = -r6/b3 - R2;
M = [R1,R0];
beta = M\rhs;
beta1 = beta(1);
beta0 = beta(2);

B2 = [beta0; beta1; 1];

% A3 = (P6 + b3*B2*P2)/P3.
term = b3*conv(B2,p2);
p6pad = poly_pad_ascending(p6,max(numel(p6),numel(term)));
termpad = poly_pad_ascending(term,numel(p6pad));
[qA,rA] = poly_divide_ascending(p6pad + termpad,p3);
A3 = poly_pad_ascending(qA,4);
rA = poly_pad_ascending(rA,3);

% A3 = alpha0 + alpha1 lambda + alpha2 lambda^2 + lambda^3.
alpha0 = A3(1);
alpha1 = A3(2);
alpha2 = A3(3);

% From A3 = (lambda-a3)B2 - b4(lambda-a5).
a3 = beta1 - alpha2;

Dpoly = conv([-a3;1],B2) - A3;
Dpoly = poly_pad_ascending(Dpoly,4);

% Dpoly should be b4*lambda - b4*a5.
b4 = Dpoly(2);
if b4 <= 1e-12
    % Roundoff or non-realizable state.  Keep positive for sqrt.
    b4 = max(abs(b4),1e-12);
end
a5 = -Dpoly(1)/b4;

% B2 = lambda^2 -(a4+a5)lambda + a4*a5 - b5.
a4 = -beta1 - a5;
b5 = a4*a5 - beta0;
if b5 <= 1e-12
    b5 = max(abs(b5),1e-12);
end

a = [a0; a1; a2; a3; a4; a5];
b = [b1; b2; b3; b4; b5];

pCheck = recurrence_poly_ascending(a,b);
pCheck = poly_pad_ascending(pCheck,7);
p6 = poly_pad_ascending(p6,7);
res = pCheck - p6;

info = struct();
info.residualInf = norm(res,inf);
info.residualRel = norm(res,2)/max(1,norm(p6,2));
info.remainderB2 = norm(M*beta-rhs,2);
info.remainderA3 = norm(rA,2);
info.D = D;
info.bmin = min(b);
info.success = isfinite(info.residualRel) && info.residualRel < 1e-8 && all(b > 0);
info.beta0 = beta0;
info.beta1 = beta1;

end
