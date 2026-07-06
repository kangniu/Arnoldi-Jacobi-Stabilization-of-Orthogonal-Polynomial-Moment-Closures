function verify_B6_formula_equilibrium_and_hyperbolicity()
%VERIFY_B6_FORMULA_EQUILIBRIUM_AND_HYPERBOLICITY
% Quick checks for the corrected B6 closing flux and characteristic speeds.

fprintf('\nB6 formula verification\n');

states = [
    0.0, 3.0, 0.0;
    0.2, 3.1, 0.1;
   -0.2, 3.2,-0.1;
    0.5, 3.6, 0.3;
   -0.5, 3.8,-0.3;
    1.0, 4.5, 0.5;
   -1.0, 4.8,-0.5
];

fprintf('%9s %9s %9s %12s %12s %12s\n', ...
    'W3','W4','W5','W6','maxImagFD','maxSpeedFD');

for i = 1:size(states,1)
    W3 = states(i,1); W4 = states(i,2); W5 = states(i,3);
    U = [1;0;1;W3;W4;W5];

    U6 = closure_raw_moment(U,'B6');
    Jfd = flux_jacobian_fd(U,'B6');
    lam = eig(Jfd);

    fprintf('%9.3f %9.3f %9.3f %12.5f %12.4e %12.4e\n', ...
        W3,W4,W5,U6,max(abs(imag(lam))),max(abs(real(lam))));
end

end
