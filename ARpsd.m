% This function returns an anonymous function that
% computes the power spectral density of an AR(p) process of order p. The
% output function will accept as input an innovations variance (S), a
% vector of lag coefficients (p), a vector of frequencies at which to
% evaluate the density (f), and the Nyquist frequency (fn)
% using modified definition from p. 392 of Percival and Walden:
% S(f) = 1/fn * sig^2/|1 - sum(p_j exp(-ij pi f/fn))|^2 

function psd = ARpsd(p)
    
defval('p',1);
% need to add recursive terms as nested functions in summed term
sumd = @(rho,f,fn) rho(1)*exp(-1i*pi*f/fn); % first term in sum
for j = 2:p
    tmp = @(rho,f,fn) rho(j)*exp(-1i*j*pi*f/fn);
    sumd = @(rho,f,fn) sumd(rho,f,fn) + tmp(rho,f,fn);
end
psd = @(S,rho,f,fn) 1/fn * S * 1./abs(1-sumd(rho,f,fn)).^2;

end