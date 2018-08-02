% This function generates confidence intervals for the power spectrum of a 
% pink noise process paramterized by an exponent specifying the power law 
% for the given process as well as a target variance of the time series 
% corresponding to the process. By default returns the 95 percentile. 
% Current version only uses pmtm to compute the power spectral density 
% estimates with a given time half-bandwidth product.
%
% IN:
% A: power law exponent for pink noise process
% var: finite variance of process
% 'conf': confidence intervals to evaluate at, can be a vector. Percents, 
%   not decimals
% 'ntrial': number of monte carlo trials (default 1000)
% 'nsample': number of data points for simulated time series (default 1000)
% 'dt': sample spacing (default 1)
% 'nw': time half-bandwidth product for multi-taper estimates (default 2)
% 'estimator': which spectral estimator to use when constructing confidence
%   intervals. Either 'pmtm' or 'pchave' (default 'pmtm'). When using
%   pchave, the 'dpss' window is chosen by default.
%
% OUT:
% CI: confidence interval(s) corresponding to requested percentiles,
%   nf x npercent
% w: frequency axis for the confidence intervals
%
% TO DO:
% - allow for various windows when using pchave
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 25.07.2018

function [CI,w] = pinkconf(A,var,varargin)

parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'A',@(a) a > 0 && a < 2);
addRequired(parser,'var',validScalarPosNum)
addParameter(parser,'nsample',1000,validScalarPosNum);
addParameter(parser,'ntrial',1000,validScalarPosNum);
addParameter(parser,'conf',95,@(x) all([x(:) > 0; x(:) < 100]));
addParameter(parser,'dt',1,@isscalar);
addParameter(parser,'nw',2,@isscalar);
addParameter(parser,'estimator','pmtm',@ischar);

parse(parser,A,var,varargin{:});
    
A      = parser.Results.A;
varnce = parser.Results.var;
n      = parser.Results.nsample;
nt     = parser.Results.ntrial;
conf   = parser.Results.conf;
dt     = parser.Results.dt;
nw     = parser.Results.nw;
est    = parser.Results.estimator;

% validate estimator
est = validatestring(est,{'pmtm','pchave'});

% compute their power spectral densities
switch est
    case 'pmtm'
        % generate t pink noise instances
        ts = pinknoise(A,n,'ntrial',nt,'ncoeff',1000,'var',varnce);
        [pxx,w] = pmtm(ts,nw,n,1/dt);
    case 'pchave'
        % generate t pink noise instances, but need to make slightly longer
        % so that we can have overlap in pchave
        ts = pinknoise(A,round(1.5*n),'ntrial',nt,'ncoeff',1000,...
            'var',varnce);
        pxx = zeros(n+1,nt);
        for ii = 1:nt
            tscell = {ts(:,ii)};
            [pxx(:,ii),w] = pchave({tscell{:}},n,95,2*n,1/dt,[],'dpss',nw);
        end
end

% get requested percentiles
CI = prctile(pxx',conf);

end