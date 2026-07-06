function nMom = closure_num_moments(closure)
switch upper(closure)
    case 'B4'
        nMom = 4;  % U0,...,U3
    case 'B6'
        nMom = 6;  % U0,...,U5
    case 'B8'
        nMom = 8;  % U0,...,U7
    case 'HYQMOM5'
        nMom = 5;  % U0,...,U4
    otherwise
        error('Unknown closure %s',closure);
end
end
