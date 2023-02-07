function w = myLog( z, tau )
% myLog: branch of the logarithm specified by
tau

% First get principal value of the log:
w = log(z);

argz = imag(w);

% Adjust if necessary...
temp = find(argz <= tau); size(temp)
while (~isempty(temp))
    argz(temp) = argz(temp)+2*pi;
    temp = find(argz <= tau);
     size(temp)
end

temp = find(argz > tau+2*pi); size(temp)
while (~isempty(temp))
    argz(temp) = argz(temp) - 2*pi;
    temp = find(argz > tau+2*pi);
     size(temp)
end

w = real(w) + 1i*argz;

end

