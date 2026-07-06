function compare_B8_static_original_vs_jacobi()
%COMPARE_B8_STATIC_ORIGINAL_VS_JACOBI Static B8 wave-speed comparison.

outdir = case_output_dir();

states = [
    0.0, 3.0, 0.0, 15.0, 0.0;
    0.2, 3.1, 0.1, 16.0, 0.8;
   -0.2, 3.2,-0.1, 17.0,-0.8;
    0.5, 3.6, 0.3, 24.0, 2.5;
   -0.5, 3.8,-0.3, 26.0,-2.5;
    1.0, 4.8, 0.5, 45.0, 6.0
];

fprintf('\nB8 static wave-speed comparison: original Jacobian vs V+A/Jacobi\n');
fprintf('%8s %8s %8s %8s %8s %12s %12s %12s %12s %9s\n', ...
    'W3','W4','W5','W6','W7','errInf','maxOrig','maxJac','resRel','bmin');

err = zeros(size(states,1),1);

for i = 1:size(states,1)
    W3 = states(i,1); W4 = states(i,2); W5 = states(i,3);
    W6 = states(i,4); W7 = states(i,5);
    U = [1;0;1;W3;W4;W5;W6;W7];

    lamOrig = sort(wave_speeds_by_jacobian(U,'B8'));
    [lamJac,info] = wave_speeds_B8_jacobi(U);
    lamJac = sort(lamJac);

    err(i) = max(abs(lamOrig-lamJac));
    fprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %12.4e %12.4e %12.4e %12.4e %9.2e\n', ...
        W3,W4,W5,W6,W7,err(i),max(abs(lamOrig)),max(abs(lamJac)), ...
        info.residualRel,info.bmin);
end

figure('Name','B8 static original vs Jacobi');
semilogy(1:numel(err),err+eps,'-o','LineWidth',1.5);
grid on;
xlabel('Test state index');
ylabel('max absolute wave-speed difference');
title('B8 static wave speeds: original Jacobian vs V+A/Jacobi');
drawnow; saveas(gcf,fullfile(outdir,'B8_static_wave_speed_original_vs_jacobi.png'));
end
