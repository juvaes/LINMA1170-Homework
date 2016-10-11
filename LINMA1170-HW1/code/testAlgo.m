function testAlgo()

clc;
clear all;
close all;

P = [1 0 -10 -7 27 70 -20 -189 50 140 0 -350];
roots(P);
alpha = 0;
beta = 1000;

[V, R] = AlgEuclideJulien(P,alpha,beta);
Multipli = ones(1,V); % vecteur contenant les multiplicités de chaque racine réelle
[coeff, degre] = firstNonZeroCoefficient(R);
while degre>0
    [V, R] = AlgEuclideJulien(R,alpha,beta)
    Multipli(1:V) = Multipli(1:V)+1;
    [coeff, degre] = firstNonZeroCoefficient(R);
end

Result = Multipli

end

% Renvoie le premier coefficient non nul d'un vecteur
function [coeff, degre] = firstNonZeroCoefficient(vec)
coeff = NaN;
degre = length(vec)-1;
for i = 1:length(vec)
    if vec(i) ~= 0
        coeff = vec(i);
        break;
    end
    degre = degre -1;
end
end



