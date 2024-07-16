function bool = isodd(number)
% ISODD  Returns true for odd and false for even input numbers
%
% 1/2010 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~isfloat(number), number = double(number); end
if number-floor(number) ~= 0
    warning('Input number not integer')
end

bool = logical(number/2-floor(number/2) ~= 0);

end      % main function isodd
