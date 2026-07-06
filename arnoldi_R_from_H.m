function R = arnoldi_R_from_H(H,n)
%ARNOLDI_R_FROM_H Build the coordinate matrix R such that V ~= Q R.
%
% If the Arnoldi relation is x*q_k = sum_j H(j,k) q_j, then the monomial
% powers can be represented as powers of the finite Hessenberg coordinate
% operator.  This routine constructs R column-by-column by reproducing
% the monomial recurrence in the Arnoldi coordinate basis.

R = zeros(n+1,n+1);
R(1,1) = 1;

for k = 1:n
    % x * previous polynomial in Arnoldi coordinates.
    % The finite recurrence is encoded by H(:,1:n).
    c = R(:,k);
    newc = zeros(n+1,1);
    for col = 1:k
        if abs(c(col)) > 0
            newc(1:col+1) = newc(1:col+1) + c(col)*H(1:col+1,col);
        end
    end
    R(:,k+1) = newc;
end

end
