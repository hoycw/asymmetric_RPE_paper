function [rho, rho_ortho] = rh_powcorrortho(A, B)

% orthogonalized power correlations as described in Hipp et al. 2012
% but based on Hilbert-transformed data and not wavelets as in the
% original version

% input are two complex-valued time series of equal length
% 1xN length

% written by
% rhelfrich@berkeley.edu
% 06-04-2018

% orthogonalized signals A' and B'
Bortho = imag(B .* (conj(A) ./ abs(A)));
Aortho = imag(A .* (conj(B) ./ abs(B)));

% log-transform of power values
Apow = log10(abs(A).^2);
Bpow = log10(abs(B).^2);
Aorthopow = log10(Aortho.^2);
Borthopow = log10(Bortho.^2);

% correlation
rho = diag(corr(Apow', Bpow'));

% correlate A' to B and B' to A;
r1 = diag(corr(Aorthopow', Bpow'));
r2 = diag(corr(Borthopow', Apow'));

% avg the correlations
rho_ortho = (r1+r2) ./ 2;

end