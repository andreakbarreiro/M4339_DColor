function w = EssSingTrunc( z,k )
% EssSingTrunc: visualie what happens near an essential singularity
%
%   exp(1/z), truncated at order k
%  

% If necessary...
k = round(k);

w = 1;

for k1=1:k
    w = w + 1./(z.^k1)/factorial(k1);
end

end

