function verify_B8_formula_equilibrium_and_hyperbolicity()
%VERIFY_B8_FORMULA_EQUILIBRIUM_AND_HYPERBOLICITY Basic B8 closure checks.

states = [
    0.0, 3.0, 0.0, 15.0, 0.0;
    0.2, 3.1, 0.1, 16.0, 0.8;
   -0.2, 3.2,-0.1, 17.0,-0.8;
    0.5, 3.6, 0.3, 24.0, 2.5;
   -0.5, 3.8,-0.3, 26.0,-2.5;
    1.0, 4.8, 0.5, 45.0, 6.0
];

fprintf('\nB8 formula verification\n');
fprintf('%8s %8s %8s %8s %8s %12s %12s %12s %12s\n', ...
    'W3','W4','W5','W6','W7','W8','maxImagFD','maxSpeedFD','polyAD-FD');

for i = 1:size(states,1)
    W3 = states(i,1); W4 = states(i,2); W5 = states(i,3);
    W6 = states(i,4); W7 = states(i,5);
    W8 = b8_closing_moment_normalized(W3,W4,W5,W6,W7);
    U = [1;0;1;W3;W4;W5;W6;W7];
    J = flux_jacobian_fd(U,'B8');
    lam = eig(J);
    pAD = b8_charpoly_coeffs_normalized(W3,W4,W5,W6,W7);
    pFD = b8_charpoly_coeffs_normalized_fd(W3,W4,W5,W6,W7);
    polyErr = norm(pAD-pFD,2)/max(1,norm(pAD,2));
    fprintf('%8.3f %8.3f %8.3f %8.3f %8.3f %12.5f %12.4e %12.4e %12.4e\n', ...
        W3,W4,W5,W6,W7,W8,max(abs(imag(lam))),max(abs(lam)),polyErr);
end
end
