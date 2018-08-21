% This function returns the coefficients for a digital filter that operates
% on white gaussian noise to generate pink noise, as described by Kasdin
% 1995. The user can choose whether to use a moving average (MA) filter,
% with coefficients computed by eq. 103-104 (pg. 821) or an autoregressive
% (AR) filter, with coefficients computed by eq. 116 (pg. 822).
% Autoregressive filter coefficients are returned by default, which can be
% fed to matlab's filter function as the second argument (the moving
% average coefficients would then be the first argument).
% The pink noise will have a power law slope of -a (so a is positive as 
% input).
%
% IN:
% A: power law coefficient ([0,2])
% 'ncoeff': (default 50) number of coefficients to generate 
% 'filter': (default 'ar') type of filter coefficients to generate, either
%   moving average 'ma' or autoregressive 'ar'
%
% OUT:
% a: filter coefficients, will be ncoeff+1 long because the first entry is
%   always equal to one
%
% TO DO:
%
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 21.08.2018

function a = pinkcoeff(A,varargin)

% parse inputs
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'A',@(A) A >= 0 && A <= 2)
addOptional(parser,'ncoeff',50,validScalarPosNum)
addParameter(parser,'filter','ar',@ischar)
parse(parser,A,varargin{:})

A = parser.Results.A;
ncoeff = parser.Results.ncoeff;
filttype = parser.Results.filter;

% validate filter
filttype = validatestring(filttype,{'ma','ar'});

% generate coefficients
a = ones(ncoeff+1,1);
switch filttype
    case 'ar'
        for ii = 1:ncoeff
            a(ii+1) = (ii - 1 - A/2)*(a(ii)/ii);
        end
    case 'ma'
        for ii = 1:ncoeff
            a(ii+1) = (A/2 + ii - 1)*(a(ii)/ii);
        end
end