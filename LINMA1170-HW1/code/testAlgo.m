function testAlgo()
clc;
clear all;
close all;

LASTN = maxNumCompThreads(1)

P = [1 0 -10 -7 27 70 -20 -189 50 140 0 -350];
r = roots(P)
alpha = -100;
beta = 100;

[V, R] = AlgEuclide(P,alpha,beta);
Multipli = ones(1,V); % vecteur contenant les multiplicités de chaque racine réelle
[coeff, degre] = firstNonZeroCoefficient(R);
while degre>0
    [V, R] = AlgEuclide(R,alpha,beta);
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
    if abs(vec(i))>10^(-6)
        coeff = vec(i);
        break;
    else
        degre = degre -1;
    end
end
end



