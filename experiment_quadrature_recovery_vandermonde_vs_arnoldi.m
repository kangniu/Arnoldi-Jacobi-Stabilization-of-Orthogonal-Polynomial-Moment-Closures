function experiment_quadrature_recovery_vandermonde_vs_arnoldi()
%EXPERIMENT_QUADRATURE_RECOVERY_VANDERMONDE_VS_ARNOLDI
% Recover positive quadrature weights using two representations.
%
% The key V+A principle is NOT to solve an ill-conditioned Vandermonde system
% and NOT to convert a high-order monomial representation at the end.  Instead
% the computation is formulated in an Arnoldi-generated basis from the
% beginning.
%
% Monomial route:
%     U = V^T w,      recover w from V^T w = U.
%
% Arnoldi-basis route:
%     qMom = Q^T w,   recover w from Q^T w = qMom.
%
% For reference we also report the unsafe post-processing route
%     qMom_bad = C^T U,
% where C contains monomial coefficients of the Arnoldi polynomials.  This
% route can be unstable because C itself may be badly scaled.  It is included
% only to demonstrate why V+A must be used as a basis formulation from the
% start.

outdir = case_output_dir();

orders = [10 20 30 40 50];
scales = [1, 1e-2, 1e-4];

rng(1);

fprintf('\nQuadrature weight recovery: monomial Vandermonde vs Arnoldi basis\n');
fprintf('%6s %10s %15s %15s %13s %13s %13s %9s\n', ...
    'N','scale','cond(V^T)','cond(Q^T)','relerr V','relerr A','A bad','neg V');

ERRV = zeros(numel(orders),numel(scales));
ERRA = zeros(numel(orders),numel(scales));
ERRBAD = zeros(numel(orders),numel(scales));

warningState = warning('off','MATLAB:nearlySingularMatrix');

for is = 1:numel(scales)
    scale = scales(is);
    for io = 1:numel(orders)
        N = orders(io);

        J = legendre_jacobi_matrix(N,1,scale);
        nodes = sort(eig((J+J')/2));

        wTrue = 0.1 + rand(N,1);
        wTrue = wTrue/sum(wTrue);

        V = zeros(N,N);
        V(:,1) = 1;
        for k = 2:N
            V(:,k) = nodes .* V(:,k-1);
        end

        % Monomial moments.
        U = V' * wTrue;

        % Original monomial Vandermonde recovery.
        wV = V' \ U;

        % Arnoldi basis built on the same nodes.
        [Q,~,C] = arnoldi_basis_values_coeffs(nodes,N-1,true);

        % Correct V+A-style representation: data are represented in the
        % Arnoldi basis directly, so no ill-conditioned monomial solve or
        % monomial-to-Arnoldi coefficient conversion is required.
        qMom = Q' * wTrue;
        wA = Q' \ qMom;

        % Unsafe post-processing: convert monomial moments to Arnoldi moments
        % by using monomial coefficients. This is intentionally diagnostic.
        qMomBad = C' * U;
        wBad = Q' \ qMomBad;

        relV = norm(wV-wTrue)/norm(wTrue);
        relA = norm(wA-wTrue)/norm(wTrue);
        relBad = norm(wBad-wTrue)/norm(wTrue);
        negV = sum(wV < -1e-10);

        ERRV(io,is) = relV;
        ERRA(io,is) = relA;
        ERRBAD(io,is) = relBad;

        fprintf('%6d %10.1e %15.4e %15.4e %13.4e %13.4e %13.4e %9d\n', ...
            N,scale,cond(V'),cond(Q'),relV,relA,relBad,negV);
    end
end

warning(warningState);

figure('Name','Quadrature recovery corrected comparison');
for is = 1:numel(scales)
    semilogy(orders,ERRV(:,is)+eps,'-o','LineWidth',1.3, ...
        'DisplayName',['V scale=',num2str(scales(is),'%.0e')]);
    hold on;
    semilogy(orders,ERRA(:,is)+eps,'--s','LineWidth',1.3, ...
        'DisplayName',['Arnoldi scale=',num2str(scales(is),'%.0e')]);
end
grid on;
xlabel('Number of nodes/moments N');
ylabel('relative weight error');
title('Weight recovery: ill-conditioned monomial system vs Arnoldi basis');
legend('Location','bestoutside');
drawnow; saveas(gcf,fullfile(outdir,'quadrature_weight_error_corrected.png'));

figure('Name','Unsafe monomial-to-Arnoldi conversion');
for is = 1:numel(scales)
    semilogy(orders,ERRBAD(:,is)+eps,'-^','LineWidth',1.3, ...
        'DisplayName',['bad transform scale=',num2str(scales(is),'%.0e')]);
    hold on;
end
grid on;
xlabel('Number of nodes/moments N');
ylabel('relative error');
title('Diagnostic: post-hoc monomial-to-Arnoldi conversion is not V+A');
legend('Location','northwest');
drawnow; saveas(gcf,fullfile(outdir,'quadrature_bad_transform_diagnostic.png'));

end
