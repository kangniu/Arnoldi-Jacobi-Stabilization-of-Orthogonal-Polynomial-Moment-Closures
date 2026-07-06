function experiment_high_order_wave_speed_stability()
%EXPERIMENT_HIGH_ORDER_WAVE_SPEED_STABILITY
% Compare high-order wave-speed computation by:
%   Original monomial route: coefficients poly(J), then roots(...)
%   Arnoldi/Jacobi route:   eig(J)
%
% We use shifted/scaled Legendre Jacobi matrices to emulate characteristic
% polynomials whose roots are real but increasingly clustered.  This isolates
% the numerical linear algebra issue that appears in high-order moment
% closures: monomial characteristic polynomials are sensitive, while symmetric
% Jacobi eigenproblems are stable.

outdir = case_output_dir();

orders = [10 20 30 40 60 80];
scales = [1, 1e-2, 1e-4];

fprintf('\nHigh-order wave-speed stability: roots(poly(J)) vs eig(J)\n');
fprintf('%6s %10s %15s %15s %15s %15s\n', ...
    'N','scale','errInf roots','errRMS roots','maxIm roots','cond(Vroots)');

ERR = zeros(numel(orders),numel(scales));

for is = 1:numel(scales)
    scale = scales(is);
    for io = 1:numel(orders)
        N = orders(io);
        J = legendre_jacobi_matrix(N,1,scale);
        lambdaJ = sort(eig((J+J')/2));

        % Original monomial/companion-style route.
        c = poly(J);
        lambdaRoots = roots(c);

        [errInf,errRMS,maxImag] = sorted_root_error(lambdaRoots,lambdaJ);

        % Vandermonde conditioning of the true roots, as a diagnostic.
        V = zeros(N,N);
        V(:,1) = 1;
        for k = 2:N
            V(:,k) = lambdaJ .* V(:,k-1);
        end
        condV = cond(V);

        ERR(io,is) = errInf;

        fprintf('%6d %10.1e %15.4e %15.4e %15.4e %15.4e\n', ...
            N,scale,errInf,errRMS,maxImag,condV);
    end
end

figure('Name','High-order wave-speed stability');
for is = 1:numel(scales)
    semilogy(orders,ERR(:,is)+eps,'-o','LineWidth',1.5, ...
        'DisplayName',['cluster scale=',num2str(scales(is),'%.0e')]);
    hold on;
end
grid on;
xlabel('Matrix/order N');
ylabel('max root error');
title('Monomial roots vs symmetric Jacobi eigenvalues');
legend('Location','northwest');
drawnow; saveas(gcf,fullfile(outdir,'high_order_roots_vs_jacobi.png'));

end
