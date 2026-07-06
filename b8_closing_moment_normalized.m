function [W8,info] = b8_closing_moment_normalized(W3,W4,W5,W6,W7)
%B8_CLOSING_MOMENT_NORMALIZED Compute the B8 normalized closing moment.
%
% The paper's B8 hierarchy sets
%   b4 = 2/3*(b1+b2+b3)
%        + 1/3*((a3-a0)^2+(a3-a1)^2+(a3-a2)^2).
% Since b4=<P4^2>/<P3^2> is affine in W8, W8 is recovered by one scalar
% solve instead of hand-transcribing a very large symbolic expression.

Wknown = [1;0;1;W3;W4;W5;W6;W7;0];

[a3v,b3v,polys,norms] = recurrence_from_moments(Wknown,4);
a0 = a3v(1); a1 = a3v(2); a2 = a3v(3); a3 = a3v(4);
b1 = b3v(1); b2 = b3v(2); b3 = b3v(3);

b4Target = (2/3)*(b1+b2+b3) + ...
    (1/3)*((a3-a0)^2 + (a3-a1)^2 + (a3-a2)^2);

p4 = polys{5};
if abs(p4(end)-1) > 1e-10
    error('P4 is expected to be monic.');
end

% <P4,P4> = W8 + terms involving W0,...,W7.
rest = moment_poly_inner(p4,p4,Wknown);
normP3 = norms(4);
W8 = b4Target*normP3 - rest;

info = struct();
info.aKnown = [a0;a1;a2;a3];
info.bKnown = [b1;b2;b3;b4Target];
info.b4Target = b4Target;
info.normP3 = normP3;
info.p4 = p4;
info.margin = min([b1,b2,b3,b4Target]);

if nargout > 1
    [gradW8,gradInfo] = b8_closing_moment_gradient(W3,W4,W5,W6,W7);
    info.gradW8 = gradW8(:);
    info.gradInfo = gradInfo;
end
end

function [gradW8,gradInfo] = b8_closing_moment_gradient(W3,W4,W5,W6,W7)
% Sensitivity of W8 with respect to [W3,W4,W5,W6,W7].

Wknown = [1;0;1;W3;W4;W5;W6;W7;0];
dW = zeros(numel(Wknown),5);
dW(4,1) = 1;
dW(5,2) = 1;
dW(6,3) = 1;
dW(7,4) = 1;
dW(8,5) = 1;

[a,da,b,db,polys,dpolys,norms,dnorms] = recurrence_from_moments_grad(Wknown,dW,4);
a0 = a(1); a1 = a(2); a2 = a(3); a3 = a(4);
b1 = b(1); b2 = b(2); b3 = b(3);

da0 = da(1,:); da1 = da(2,:); da2 = da(3,:); da3 = da(4,:);
db1 = db(1,:); db2 = db(2,:); db3 = db(3,:);

b4Target = (2/3)*(b1+b2+b3) + ...
    (1/3)*((a3-a0)^2 + (a3-a1)^2 + (a3-a2)^2);

db4Target = (2/3)*(db1+db2+db3) + ...
    (2/3)*(a3-a0)*(da3-da0) + ...
    (2/3)*(a3-a1)*(da3-da1) + ...
    (2/3)*(a3-a2)*(da3-da2);

p4 = polys{5};
dp4 = dpolys{5};
[rest,drest] = moment_poly_inner_grad(p4,dp4,p4,dp4,Wknown,dW);
normP3 = norms(4);
dnormP3 = dnorms(4,:);

gradW8 = db4Target*normP3 + b4Target*dnormP3 - drest;

gradInfo = struct();
gradInfo.b4Target = b4Target;
gradInfo.db4Target = db4Target(:);
gradInfo.rest = rest;
gradInfo.drest = drest(:);
gradInfo.normP3 = normP3;
gradInfo.dnormP3 = dnormP3(:);
end
