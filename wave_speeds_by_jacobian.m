function lambda = wave_speeds_by_jacobian(U,closure)
%WAVE_SPEEDS_BY_JACOBIAN Finite-difference flux Jacobian eigenvalues.
J = flux_jacobian_fd(U,closure);
lambda = eig(J);
lambda = sort(real(lambda));
end
