% This function returns an anonymous function that computes the power
% spectral density of a pink noise process. The output function will accept
% as input a power law exponent A, power law coefficient C, and a vector of
% frequencies at which to evaluate the density (f).

function psd = pinkpsd()

psd = @(A,C,f) C*f.^(-A);
    
end