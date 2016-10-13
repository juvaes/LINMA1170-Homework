close all; clc;

alpha = -2;
beta = 2;

e = 10^(-41);
coeff = [1 (1-e) e 2*e*(1-e) -2*e*e];

y = polyval(coeff,-1)


[V,R] = AlgEuclideJulien(coeff,alpha,beta);
s = V
r = roots(coeff)
rend = r(end)

