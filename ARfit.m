% This function fits the innovations variance and lag coefficients of an 
% AR(p) process to the observed power spectral density estimate of a time
% series.
%
% IN:
% p: order of AR process to consider
% f: frequency axis of time series
% pxx: power spectral density estimate
% fn: Nyquist frequency of time series
% 'S0' (optional): estimate of innovations variance, estimated as total
%   signal variance by default
% 'rho0' (optional): estimate of lag coefficient(s), computed by invfreqz
%   by default
% 'weights': (default 1/sqrt(f)) apply weighting when computing the
%   best-fitting AR(p) spectrum. weights must be of equal length to f, pxx
%
% OUT:
% rho: optimal lag coefficient(s)
% S: optimal innovations variance
%
% TO DO:
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 20.08.2018

function [rho,S] = ARfit(p,f,pxx,fn,varargin)
    
parser = inputParser;
addRequired(parser,'p',@isscalar)
addRequired(parser,'f',@isnumeric)
addRequired(parser,'pxx',@isnumeric)
addRequired(parser,'fn',@isscalar)
addParameter(parser,'S0',[],@isnumeric)
addParameter(parser,'rho0',[],@isnumeric)
addParameter(parser,'weights',1./sqrt(f),@isnumeric)

% parse inputs
parse(parser,p,f,pxx,fn,varargin{:});
p     = parser.Results.p;
f     = parser.Results.f;
pxx   = parser.Results.pxx;
fn    = parser.Results.fn;
S0    = parser.Results.S0;
rho0  = parser.Results.rho0;
wghts = parser.Results.weights;

% make sure that given frequencies are of same length as psd of data
assert(length(f) == length(pxx),'f and pxx must be same length')

% make sure weights are same length as psd of data
assert(length(wghts) == length(pxx),'weights and pxx must be same length')

% generate initial guess for rho for fminsearch
if isempty(rho0)
    [~,rho0] = invfreqz(sqrt(pxx),linspace(0,pi,length(f)),0,p,[],100);
    rho0 = -rho0(2:end);
else
    assert(length(rho0) == p,'rho0 must be of length equal to p')
end

% generate initial guess for S0 for fminsearch
if isempty(S0)
    S0 = sum(pxx)*mean(diff(f)); % make variance of process
else
    assert(length(S0) == 1,'S0 must be a single value')
end

% ensure rho0 is column
rho0 = rho0(:);

% ignore zero frequency
idx = f~=0;
f = f(idx);
pxx = pxx(idx);

% get functional form of AR(p) power spectral density
psd = ARpsd(p);
% generate objective function
obj = @(x0) sum( wghts.*abs( log(psd(x0(1),x0(2:end),f,fn))- log(pxx) ) );
% look for optimal S, rho starting at S0, rho0
%     opt = optimset('MaxFunEvals',600*p,'MaxIter',600*p);
nonlcon = @lagroots; 
X = fmincon(obj,[S0;rho0],[],[],[],[],...
    [],[],nonlcon);
S = X(1);
rho = X(2:end); 
    
%% nonlinear constraint function for fmincon
% from pg 392, Percival and Walden. I impose the constraint that the roots 
% lie outside the unit circle, so that subtracting them from 1 imposes 
% non-positivity, which is what matlab's nonlincon can achieve
function [c,ceq] = lagroots(x)
    c = 1 - abs(roots([flipud(-x(2:end));1]));
    ceq = [];
end

end