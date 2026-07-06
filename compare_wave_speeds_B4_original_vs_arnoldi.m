function compare_wave_speeds_B4_original_vs_arnoldi()
%COMPARE_WAVE_SPEEDS_B4_ORIGINAL_VS_ARNOLDI
% Compare the original flux-Jacobian eigenvalue computation with the
% Arnoldi/Jacobi wave-speed computation for the B4 closure.

outdir = case_output_dir();

W3grid = linspace(-6,6,401);
err = zeros(size(W3grid));
imagOrig = zeros(size(W3grid));
maxSpeedOrig = zeros(size(W3grid));
maxSpeedArn = zeros(size(W3grid));

for i = 1:numel(W3grid)
    W3 = W3grid(i);

    % B4 normalized state: rho=1,u=0,theta=1, U3=W3.
    U = [1;0;1;W3];

    Jfd = flux_jacobian_fd(U,'B4');
    lamOrigFull = eig(Jfd);
    lamOrig = sort(real(lamOrigFull));
    lamArn = sort(wave_speeds_B4_jacobi(U));

    err(i) = max(abs(lamOrig-lamArn));
    imagOrig(i) = max(abs(imag(lamOrigFull)));
    maxSpeedOrig(i) = max(abs(lamOrig));
    maxSpeedArn(i) = max(abs(lamArn));
end

fprintf('\nB4 wave-speed comparison: original Jacobian vs Arnoldi/Jacobi\n');
fprintf('  max |lambda_original - lambda_arnoldi| = %.4e\n',max(err));
fprintf('  max imaginary part of original Jacobian eigs = %.4e\n',max(imagOrig));
fprintf('  max speed original = %.4e, Arnoldi = %.4e\n',max(maxSpeedOrig),max(maxSpeedArn));

figure('Name','B4 wave-speed comparison');
semilogy(W3grid,err + eps,'LineWidth',1.5); grid on;
xlabel('W_3');
ylabel('max absolute eigenvalue difference');
title('B4: flux-Jacobian eigenvalues vs Arnoldi/Jacobi eigenvalues');
drawnow; saveas(gcf,fullfile(outdir,'B4_wave_speed_error_original_vs_arnoldi.png'));

figure('Name','B4 original imaginary parts');
semilogy(W3grid,imagOrig + eps,'LineWidth',1.5); grid on;
xlabel('W_3');
ylabel('max imaginary part');
title('B4: imaginary part from numerical flux-Jacobian eigenvalues');
drawnow; saveas(gcf,fullfile(outdir,'B4_original_jacobian_imaginary_parts.png'));

end
