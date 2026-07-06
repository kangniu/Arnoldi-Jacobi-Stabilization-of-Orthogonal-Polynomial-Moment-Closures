function [errInf,errRMS,maxImag] = sorted_root_error(z,lambdaTrue)
%SORTED_ROOT_ERROR Compare computed roots with sorted true real eigenvalues.
%
% This assumes a one-dimensional ordered spectrum.  The real parts are sorted
% before comparison.

z = z(:);
lambdaTrue = lambdaTrue(:);
zr = sort(real(z));
lt = sort(real(lambdaTrue));
n = min(numel(zr),numel(lt));

err = zr(1:n) - lt(1:n);
errInf = max(abs(err));
errRMS = sqrt(mean(abs(err).^2));
maxImag = max(abs(imag(z)));
end
