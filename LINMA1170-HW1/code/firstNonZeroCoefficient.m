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