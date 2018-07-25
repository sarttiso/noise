% This function returns the coefficients for a moving average filter that
% operates on white gaussian noise to generate pink noise, as described by
% Kasdin 1995, pg 821, eq 103. the pink noise will have a power law slope
% of -a (so a is positive as input)
%
% IN:
% ncoeff: number of coefficients to generate (default 50)
%
% OUT:
% b: filter coefficients

function b = pinkcoeff(a,varargin)

parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'a',validScalarPosNum)
addOptional(parser,'ncoeff',50,validScalarPosNum)
parser(parser,a,varargin{:})

a = parser.Results.a;
ncoeff = parser.Results.ncoeff;

b = ones(ncoeff,1);
for ii = 1:ncoeff-1
    b(ii+1) = (a/2 + ii - 1)*(b(ii)/ii);
end

end