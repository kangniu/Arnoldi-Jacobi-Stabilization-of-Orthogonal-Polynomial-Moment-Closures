function T = compare_B8_charpoly_ad_vs_fd(nRepeat)
%COMPARE_B8_CHARPOLY_AD_VS_FD Compare B8 P8 construction routes.
%
% The production route uses forward sensitivities through the recurrence.
% The finite-difference route is retained only as a diagnostic reference.

if nargin < 1 || isempty(nRepeat)
    nRepeat = 50;
end

states = [
    0.0, 3.0, 0.0, 15.0, 0.0;
    0.2, 3.1, 0.1, 16.0, 0.8;
   -0.2, 3.2,-0.1, 17.0,-0.8;
    0.5, 3.6, 0.3, 24.0, 2.5;
   -0.5, 3.8,-0.3, 26.0,-2.5;
    1.0, 4.8, 0.5, 45.0, 6.0
];

rows = cell(size(states,1),1);

warm = states(1,:);
b8_charpoly_coeffs_normalized(warm(1),warm(2),warm(3),warm(4),warm(5));
b8_charpoly_coeffs_normalized_fd(warm(1),warm(2),warm(3),warm(4),warm(5));

fprintf('\nB8 characteristic polynomial: forward sensitivity vs finite difference\n');
fprintf('%4s %10s %10s %12s %12s %10s\n', ...
    'idx','errInf','errRel','timeAD','timeFD','speedup');

for i = 1:size(states,1)
    s = states(i,:);

    pAD = b8_charpoly_coeffs_normalized(s(1),s(2),s(3),s(4),s(5));
    pFD = b8_charpoly_coeffs_normalized_fd(s(1),s(2),s(3),s(4),s(5));
    errInf = norm(pAD-pFD,inf);
    errRel = norm(pAD-pFD,2)/max(1,norm(pAD,2));

    tic;
    for k = 1:nRepeat
        b8_charpoly_coeffs_normalized(s(1),s(2),s(3),s(4),s(5));
    end
    timeAD = toc/nRepeat;

    tic;
    for k = 1:nRepeat
        b8_charpoly_coeffs_normalized_fd(s(1),s(2),s(3),s(4),s(5));
    end
    timeFD = toc/nRepeat;
    speedup = timeFD/timeAD;

    fprintf('%4d %10.3e %10.3e %12.4e %12.4e %10.2f\n', ...
        i,errInf,errRel,timeAD,timeFD,speedup);

    rows{i} = {i,s(1),s(2),s(3),s(4),s(5),errInf,errRel,timeAD,timeFD,speedup};
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'idx','W3','W4','W5','W6','W7','errInf','errRel','timeAD','timeFD','speedup'});

outdir = case_output_dir();
writetable(T,fullfile(outdir,'B8_charpoly_ad_vs_fd.csv'));
end
