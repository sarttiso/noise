% This function generates confidence intervals for a pink noise process
% paramterized by a coefficient and exponent specifying the power law for
% the given process. By default returns the 95 percentile. Current version
% only uses pmtm to compute the power spectral density estimates with a
% given time half-bandwidth product.
%
% IN:
% A: power law exponent for pink noise process
% C: power law coefficient for pink noise process
% 'conf': confidence intervals to evaluate at, can be a vector. Percents, 
%   not decimals
% 't': number of monte carlo trials
% 'n': number of data points to output. 
% 'dt': sample spacing
% 'nw': time half-bandwidth product for multi-taper estimates
% 'var': variance for signals to mimic
%
% OUT:
% CI: confidence interval(s) corresponding to requested percentiles,
%   nf x npercent
% w: frequency axis for the confidence intervals
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 17.07.2018

function [CI,w] = pinkconf(A,C,varargin)

parser = inputParser;
addRequired(parser,'A',@isscalar);
addRequired(parser,'C',@isscalar);
addParameter(parser,'n',1000,@isscalar);
addParameter(parser,'t',1000,@isscalar);
addParameter(parser,'conf',95,@isnumeric);
addParameter(parser,'dt',1,@isscalar);
addParameter(parser,'nw',2,@isscalar);
% addParameter(parser,'estimator','pmtm',@ischar);

parse(parser,A,C,varargin{:});
    
A = parser.Results.A;
C   = parser.Results.C;
n   = parser.Results.n;
t   = parser.Results.t;
conf = parser.Results.conf;
dt  = parser.Results.dt;
nw = parser.Results.nw;

% generate t pink noise instances


% compute their power spectral densities
[pxx,w] = pmtm(O,nw,n,1/dt);

% get requested percentiles
CI = prctile(pxx',conf);

end