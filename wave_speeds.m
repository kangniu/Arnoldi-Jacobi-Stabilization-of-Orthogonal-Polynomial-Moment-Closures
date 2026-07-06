function lambda = wave_speeds(U,closure,waveMethod)
%WAVE_SPEEDS Compute characteristic speeds.
%
% For B4, waveMethod='jacobi' uses the Jacobi matrix implied by the
% orthogonal-polynomial recurrence.  Otherwise a finite-difference flux
% Jacobian is used.

if nargin < 3
    waveMethod = 'jacobian';
end

switch lower(waveMethod)
    case 'jacobi'
        if strcmpi(closure,'B4')
            lambda = wave_speeds_B4_jacobi(U);
        elseif strcmpi(closure,'B6')
            lambda = wave_speeds_B6_jacobi(U);
        elseif strcmpi(closure,'B8')
            lambda = wave_speeds_B8_jacobi(U);
        else
            lambda = wave_speeds_by_jacobian(U,closure);
        end
    case 'jacobian'
        lambda = wave_speeds_by_jacobian(U,closure);
    otherwise
        error('Unknown wave speed method %s',waveMethod);
end

end
