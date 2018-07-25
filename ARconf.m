% Monte Carlo simulation of AR(p) spectra with parameters fit to time
% series. Make sure that if the time series is not unit variance to
% specify the variance, or else the confidence intervals will not have the 
% appropriate magnitude.
%
% rho: lag autocorrelation coefficient(s), number of coefficients is order
%   of autoregressive process
% e: innovations variance
% 'conf': confidence intervals to evaluate at, can be a vector. Percents, 
%   not decimals
% 't': number of monte carlo trials
% 'n': number of data points to output. 
% 'dt': sample spacing
% 'estimator': 'pmtm' (default) or 'pchave'
function [C,w] = ARconf(rho,e,varargin)

parser = inputParser;
addRequired(parser,'rho',@isnumeric);
addRequired(parser,'e',@isscalar);
addParameter(parser,'n',1000,@isscalar);
addParameter(parser,'t',1000,@isscalar);
addParameter(parser,'conf',95,@isnumeric);
addParameter(parser,'dt',1,@isscalar);
addParameter(parser,'nw',2,@isscalar);
addParameter(parser,'estimator','pmtm',@ischar);

parse(parser,rho,e,varargin{:});

rho = parser.Results.rho;
e   = parser.Results.e;
n   = parser.Results.n;
t   = parser.Results.t;
conf = parser.Results.conf;
dt  = parser.Results.dt;
nw = parser.Results.nw;
est = parser.Results.estimator;

% validate estimator
est = validatestring(est,{'pmtm','pchave'});

% create AR(p) model
arm = arima('Constant',0,'AR',rho,'Variance',e);

% noise iterations
O = simulate(arm,n,'NumPaths',t);

% PSDs of noise
%     [pxx,w] = periodogram(O,[],[],1/dt);
if strcmp(est,'pmtm')
    [pxx,w] = pmtm(O,nw,n,1/dt);
elseif strcmp(est,'pchave')
    O = num2cell(O,1);
    pxx = zeros(floor(n/2)+1,t);
    for j = 1:t
        [pxx(:,j),w] = pchave(O{j},floor(n/4),80,n,1/dt,[],'dpss',nw);
    end
end
%     C = mean(pxx,2);
% get interval
C = prctile(pxx',conf);
%     C = C';

% reincorporate variance of data
%     C = C*S;

end