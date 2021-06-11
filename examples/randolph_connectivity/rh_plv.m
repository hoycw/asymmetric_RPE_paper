function [plv, iplv] = rh_plv(A, B)

% PLV and imagPLV based on the hilbert-transformed complex data

% input are two complex-valued time series of equal length
% 1xN length

% written by
% rhelfrich@berkeley.edu
% 08-16-2018

% calculate the plv and iplv
phaseA = angle(A);
phaseB = angle(B); 

plv = abs(mean(exp(1i*(phaseA - phaseB))));

% take the abs value of the imag value, 
% since the direction here is meaningless
iplv = abs(imag(mean(exp(1i*(phaseA - phaseB)))));


end