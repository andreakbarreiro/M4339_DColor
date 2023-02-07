function w = PolyTrunc( z,k )
% PolyTrunc: visualize what happens at the edge of the ROC
%
%   1/(1-z), truncated at order k
%  

% If necessary...
k = round(k);

w = 1;

for k1=1:k
    w = w + (z.^k1);
end

end

