function wv = velocity_quadrature_weights(v)
%VELOCITY_QUADRATURE_WEIGHTS Trapezoidal weights for a uniform velocity grid.
v = v(:);
Nv = numel(v);
if Nv < 2
    wv = 1;
    return;
end
dv = (v(end)-v(1))/(Nv-1);
wv = dv*ones(Nv,1);
wv(1) = 0.5*dv;
wv(end) = 0.5*dv;
end
