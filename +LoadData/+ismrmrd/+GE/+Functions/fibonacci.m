function outp=fibonacci(inp)
%FIBONACCI Calculate Fibonacci numbers
% 10/2014 Rolf Schulte
if (nargin<1), help(mfilename); return; end

outp = zeros(size(inp(:)));
for l1=1:length(inp(:))

    n = inp(l1);
    f1 = 1;
    f2 = 1;
    
    if n<0, error('only implemented for positive numbers'); end
    for l=3:n
        f = f1+f2;
        f1 = f2;
        f2 = f;
    end
    switch n
        case 0, f = 0;
        case 1, f = 1;
        case 2, f = 1;
    end
    outp(l1) = f;
end

outp = reshape(outp,size(inp));

end      % fibonacci.m
