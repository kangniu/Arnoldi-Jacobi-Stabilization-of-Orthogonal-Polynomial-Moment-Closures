function compare_B6_static_original_vs_jacobi()
%COMPARE_B6_STATIC_ORIGINAL_VS_JACOBI
% Static comparison of B6 wave speeds:
%   original route: finite-difference flux Jacobian,
%   prototype Arnoldi/Jacobi route: coefficient-matched Jacobi matrix.
%
% This test should be run before attempting full B6 PDE simulations, because
% B6 coefficient matching is currently a prototype and is more expensive.

outdir = case_output_dir();

states = [
    0.0, 3.0, 0.0;
    0.2, 3.1, 0.1;
   -0.2, 3.2,-0.1;
    0.5, 3.6, 0.3;
   -0.5, 3.8,-0.3;
    1.0, 4.5, 0.5;
   -1.0, 4.8,-0.5
];

fprintf('\nB6 static wave-speed comparison: original Jacobian vs prototype Jacobi\n');
fprintf('%9s %9s %9s %14s %14s %14s %14s %9s\n', ...
    'W3','W4','W5','errInf','maxOrig','maxJac','resRel','bmin');

err = zeros(size(states,1),1);

for i = 1:size(states,1)
    W3 = states(i,1); W4 = states(i,2); W5 = states(i,3);

    % Normalized raw state: rho=1, u=0, theta=1.
    U = [1;0;1;W3;W4;W5];

    lamOrig = sort(wave_speeds_by_jacobian(U,'B6'));

    [lamJac,info] = wave_speeds_B6_jacobi(U);
    lamJac = sort(lamJac);

    err(i) = max(abs(lamOrig-lamJac));

    fprintf('%9.3f %9.3f %9.3f %14.4e %14.4e %14.4e %14.4e %9.2e\n', ...
        W3,W4,W5,err(i),max(abs(lamOrig)),max(abs(lamJac)), ...
        info.residualRel,info.bmin);
end

figure('Name','B6 static original vs Jacobi');
semilogy(1:numel(err),err+eps,'-o','LineWidth',1.5);
grid on;
xlabel('Test state index');
ylabel('max absolute wave-speed difference');
title('B6 static wave speeds: original Jacobian vs prototype Jacobi');
drawnow; saveas(gcf,fullfile(outdir,'B6_static_wave_speed_original_vs_jacobi.png'));

end
