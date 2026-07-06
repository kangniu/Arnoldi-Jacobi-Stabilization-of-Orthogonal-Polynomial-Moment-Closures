function F = flux_moment(U,closure)
%FLUX_MOMENT Moment flux F=[U1,...,Un,U_{n+1}] for U=[U0,...,Un].
nMom = numel(U);
F = zeros(nMom,1);
F(1:nMom-1) = U(2:nMom);
F(nMom) = closure_raw_moment(U,closure);
end
