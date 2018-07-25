% This function returns the filter coefficients that generate a stream of
% pink noise from an input stream of white noise: Y = H*X with X being
% white noise and H being the transfer function for the pink noise filter.
% The coefficients are generated for a given pink noise process with
% coefficient C and exponent A.
% This process depends on number of filter coefficients included (deg)
%
% IN:
% A: pink noise exponent
% C: pink noise coefficient
% f: frequencies for which to find best digital filter coefficients
% 'deg': number of filter coefficients to include (must be 2 or larger),
%   default 3
%
% OUT:
% b: numerator digital filter coefficients (deg x 1)
% a: denominator digital filter coefficients (deg x 1) (first is always 1)

function [b,a] = getpinkfilt(A,C,f,varargin)

parser = inputParser;
addRequired(parser,'A',@isscalar);
addRequired(parser,'C',@isscalar);
addRequired(parser,'f',@isnumeric);
addParameter(parser,'deg',3,@(x) isscalar(x) & x >= 2)

% parse inputs
parse(parser,A,C,f,varargin{:});
A = parser.Results.A;
C = parser.Results.C;
f = parser.Results.f;
deg = parser.Results.deg;

% check frequencies
assert(all(f>=0),'frequencies must be non-negative')

% remove zero frequency
f = f(f>0);

% assume white noise with this constant value for all frequencies
const = 1;

% get functional form of pink noise power spectral density
psd = pinkpsd();

% generate objective function
% the first deg entries of x0 are a. The last deg-1 entries of x0 are b. x0
% is 2*deg-1 x 1.
obj = @(x0) sum( abs( psd(A,C,f)/const -  df(x0)) );
% formulate constraints
% lb = zeros(2*deg-1,1); % lower bounds
% ub = 10;   % upper bound (on A only)
% look for optimal A,C starting at A0,C0
b0 = rand(deg,1);
a0 = rand(deg-1,1);
X = fmincon(obj,[b0;a0]);
b = X(1:deg);
a = X(deg+1:end);
a = [1;a];

% this function computes the value of the digital filter with coefficients
% given by x0 at the frequencies f (accounting for the fact that the
% digital filter is usually defined from the z-transform, so here I have to
% be sure to distinguish between z and f: z = exp(2*pi*i*f) )
function y = df(x0)
    num = zeros(length(f),1);
    den = zeros(length(f),1);
    bcoeff = x0(1:deg);
    acoeff = [1;x0(deg+1:end)];
    z = exp(1i*2*pi*f);
    for ii = 1:deg
        num = num + bcoeff(ii)*z.^(-(ii-1));
        den = den + acoeff(ii)*z.^(-(ii-1));
    end
    y = abs(num./den);
end

end