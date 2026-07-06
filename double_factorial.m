function y = double_factorial(n)
%DOUBLE_FACTORIAL n!! with convention (-1)!!=1 and 0!!=1.
if n <= 0
    y = 1;
else
    y = prod(n:-2:1);
end
end
