clc;
close all;

t = P
P = [1 0 -10 -7 27 70 -20 -189 50 140 0 -350];
r = roots(P)
alpha = -1000;
beta = 1000;

[V, R] = AlgEuclideJulien(P,alpha,beta)
Multipli = ones(1,V); % vecteur contenant les multiplicités de chaque racine réelle
[coeff, degre] = firstNonZeroCoefficient(R);
while degre>0
    [V, R] = AlgEuclideJulien(R,alpha,beta);
    Multipli(1:V) = Multipli(1:V)+1;
    [coeff, degre] = firstNonZeroCoefficient(R);
end

Result = Multipli





