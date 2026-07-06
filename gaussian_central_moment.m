function w = gaussian_central_moment(k)
%GAUSSIAN_CENTRAL_MOMENT Standard normal central moment E[Z^k].
if mod(k,2)==1
    w = 0;
else
    w = double_factorial(k-1);
end
end
