% This function simulates pink noise by applying a moving average digital
% filter to white gaussian noise of length n with a filter containing
% ncoeff coefficients. The pink noise will have a spectrum that decays as 
% a.
% IN:
% a: exponent of pink noise, goes from 0 to 2
% n: number of samples for simulated time series
% 'nt': number of time series to simulate (default 1)
% 'ncoeff': number of coefficients in filter for time series (default 50)
% 'var': variance for the time series to achieve, (default 1)
%
% OUT:
% ts: pink noise time series
%
% TO DO:
% add some sort of variance control

function ts = pinknoise(a,n,varargin)

parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'a',@(a) a > 0 && a < 2)
addRequired(parser,'n',validScalarPosNum)
addParameter(parser,'nt',1,validScalarPosNum)
addParameter(parser,'ncoeff',50,validScalarPosNum)
addParameter(parser,'var',1,validScalarPosNum)
parse(parser,a,n,varargin{:})

a = parser.Results.a;
n = parser.Results.n;
nt = parser.Results.nt;
ncoeff = parser.Results.ncoeff;
varnce = parser.Results.var;

% generate white noise
x = wgn(n+10*ncoeff,nt,1);
% get coefficients
b = pinkcoeff(a,ncoeff);
% filter white noise
ts = filter(b,1,x);
ts = ts(end-n+1:end,:);


end