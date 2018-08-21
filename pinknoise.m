% This function simulates pink noise by applying a moving average digital
% filter to white gaussian noise of length n with a filter containing
% ncoeff coefficients. The pink noise will have a spectrum that decays as 
% a.
% IN:
% a: exponent of pink noise, goes from 0 to 2
% nsample: number of samples for simulated time series
% 'ntrial': number of time series to simulate (default 1)
% 'ncoeff': number of coefficients in filter for time series (default 100)
% 'var': variance for the time series to achieve, (default 1)
%
% OUT:
% ts: pink noise time series
%
% 
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 25.07.2018

function ts = pinknoise(A,nsample,varargin)

parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'A',@(A) A >= 0 && A <= 2)
addRequired(parser,'nsample',validScalarPosNum)
addParameter(parser,'ntrial',1,validScalarPosNum)
addParameter(parser,'ncoeff',100,validScalarPosNum)
addParameter(parser,'var',1,validScalarPosNum)
parse(parser,A,nsample,varargin{:})

A = parser.Results.A;
nsample = parser.Results.nsample;
nt = parser.Results.ntrial;
ncoeff = parser.Results.ncoeff;
varnce = parser.Results.var;

% generate white noise
x = randn(nsample+10*ncoeff,nt);
% get coefficients, use moving average filter
b = pinkcoeff(A,ncoeff,'filter','ma');
% filter white noise
ts = filter(b,1,x);
ts = ts(end-nsample+1:end,:);
% impose given variance on time series
ts = bsxfun(@minus,ts,mean(ts));
ts = bsxfun(@times,ts,sqrt(varnce./var(ts)));

end