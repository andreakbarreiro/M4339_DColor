function w = ExpTrunc( z,k )
% ExpSingTrunc: visualize series approximations to the exponential function
%
%   exp(z), truncated at order k
%  

% If necessary...
k = round(k);

w = 1;

for k1=1:k
    w = w + (z.^k1)/factorial(k1);
end

end

