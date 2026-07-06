function [val,grad] = moment_poly_inner_grad(p,dp,q,dq,W,dW)
%MOMENT_POLY_INNER_GRAD Moment inner product and forward sensitivities.
%
% p and q are ascending polynomial coefficient vectors.  dp and dq contain
% coefficient sensitivities, one variable per column.  W(k+1) is moment W_k,
% and dW(k+1,:) contains its sensitivities.

p = p(:);
q = q(:);
W = W(:);

nvar = size(dW,2);
if isempty(dp)
    dp = zeros(numel(p),nvar);
end
if isempty(dq)
    dq = zeros(numel(q),nvar);
end

val = 0;
grad = zeros(1,nvar);

for i = 1:numel(p)
    for j = 1:numel(q)
        idx = (i-1) + (j-1) + 1;
        if idx > numel(W)
            error('Need moment W_%d for polynomial inner product.',idx-1);
        end
        val = val + p(i)*q(j)*W(idx);
        grad = grad + dp(i,:)*q(j)*W(idx) + p(i)*dq(j,:)*W(idx) + ...
            p(i)*q(j)*dW(idx,:);
    end
end
end
