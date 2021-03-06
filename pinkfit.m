% This function fits the exponent and coefficient of a pink noise process
% to the observed power spectral density estimate of a time series. The
% pink noise process is defined by a power spectral density of 
% S(f) = C*f^(-A) for 0<A<2.
% 
% IN:
% f: frequency axis of time series
% pxx: power spectral density estimate
% A0 (optional): estimate of noise exponent
% C0 (optional): estimate of power law coefficient
% 'weights': (default 1/sqrt(f)) apply weighting when computing the
%   best-fitting AR(p) spectrum. weights must be of equal length to f, pxx
%
% OUT:
% A: optimal pink noise exponent
% C: optimal power law coefficient
%
% TO DO:
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 20.08.2018

function [A,C] = pinkfit(f,pxx,varargin)

parser = inputParser;
addRequired(parser,'f',@isnumeric)
addRequired(parser,'pxx',@isnumeric)
addParameter(parser,'A0',1,@isscalar)
addParameter(parser,'C0',1,@isscalar)
addParameter(parser,'weights',1./sqrt(f),@isnumeric)

% parse inputs
parse(parser,f,pxx,varargin{:});
f   = parser.Results.f;
pxx = parser.Results.pxx;
A0  = parser.Results.A0;
C0 = parser.Results.C0;
wghts = parser.Results.weights;

% make columns
pxx = pxx(:);
f = f(:);

% ignore zero frequency
idx = f~=0;
f = f(idx);
pxx = pxx(idx);
wghts = wghts(idx);

% make sure that given frequencies are of same length as psd of data
assert(length(f) == length(pxx),'f and pxx must be same length')

% get functional form of pink noise power spectral density
psd = pinkpsd();
% generate objective function
obj = @(x0) sum( wghts.*abs( log(psd(x0(1),x0(2),f)) - log(pxx) ) );
% formulate constraints: A <= 2, A >= 0, C >= 0
lb = [0;0]; % lower bounds
ub = [2,Inf];   % upper bound (on A only)
% options
options = optimoptions('fmincon','Display','none');
% look for optimal A,C starting at A0,C0
X = fmincon(obj,[A0;C0],[],[],[],[],lb,ub,[],options);
A = X(1);
C = X(2:end);

end