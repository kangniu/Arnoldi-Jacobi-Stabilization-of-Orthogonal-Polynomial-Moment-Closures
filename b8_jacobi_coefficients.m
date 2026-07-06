function [a,b,info] = b8_jacobi_coefficients(W3,W4,W5,W6,W7)
%B8_JACOBI_COEFFICIENTS Semi-explicit B8 Jacobi recurrence recovery.

[W8,closeInfo] = b8_closing_moment_normalized(W3,W4,W5,W6,W7);
W = [1;0;1;W3;W4;W5;W6;W7;W8];

[~,~,polys] = recurrence_from_moments(W,4);
aHead = closeInfo.aKnown;
bHead = closeInfo.bKnown;

p4 = polys{5};
p3 = polys{4};
p8 = b8_charpoly_coeffs_normalized(W3,W4,W5,W6,W7);

[aTail,bTail,tailInfo] = recover_jacobi_tail_from_split(p8,p4,p3,bHead(4),4);

a = [aHead; aTail];
b = [bHead; bTail];

pCheck = recurrence_poly_ascending(a,b);
pCheck = poly_pad_ascending(pCheck,9);
p8 = poly_pad_ascending(p8,9);
res = pCheck - p8;

info = struct();
info.W8 = W8;
info.closeInfo = closeInfo;
info.tailInfo = tailInfo;
info.residualInf = norm(res,inf);
info.residualRel = norm(res,2)/max(1,norm(p8,2));
info.bmin = min(b);
info.success = isfinite(info.residualRel) && info.residualRel < 1e-6 && all(b > 0);
end
