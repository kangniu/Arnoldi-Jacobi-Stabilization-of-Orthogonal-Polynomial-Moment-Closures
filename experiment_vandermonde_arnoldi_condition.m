function experiment_vandermonde_arnoldi_condition()
% Compare conditioning of monomial Vandermonde and Arnoldi basis.
%
% This reproduces the spirit of the condition-number experiment in the
% uploaded Vandermonde-with-Arnoldi draft.  We use Chebyshev--Lobatto nodes
% x_j = cos(pi*j/n), j=0,...,n.

outdir = case_output_dir();

nList = [10 20 30 40 50];
condV = zeros(size(nList));
condQ = zeros(size(nList));
condR = zeros(size(nList));
orthErr = zeros(size(nList));

fprintf('%6s %15s %15s %15s %15s\n', ...
    'n','cond(V)','cond(Q)','cond(R)','orthErr');

for ii = 1:numel(nList)
    n = nList(ii);
    x = cos(pi*(0:n)'/n);

    V = zeros(n+1,n+1);
    V(:,1) = 1;
    for k = 2:n+1
        V(:,k) = x .* V(:,k-1);
    end

    [Q,H] = arnoldi_basis_values(x,n,true);
    R = arnoldi_R_from_H(H,n);

    condV(ii) = cond(V);
    condQ(ii) = cond(Q);
    condR(ii) = cond(R);
    orthErr(ii) = norm(Q'*Q/(n+1) - eye(n+1),2);

    fprintf('%6d %15.4e %15.4e %15.4e %15.4e\n', ...
        n, condV(ii), condQ(ii), condR(ii), orthErr(ii));
end

figure('Name','Condition comparison');
semilogy(nList,condV,'-o','LineWidth',1.5); hold on;
semilogy(nList,condQ,'-s','LineWidth',1.5);
semilogy(nList,condR,'-^','LineWidth',1.5);
grid on;
xlabel('Polynomial degree n');
ylabel('2-norm condition number');
legend('Vandermonde V','Arnoldi basis Q','Arnoldi coordinate R', ...
    'Location','northwest');
title('Conditioning: monomial Vandermonde vs Arnoldi basis');
drawnow;
saveas(gcf,fullfile(outdir,'condition_vandermonde_arnoldi.png'));

end
