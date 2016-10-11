function [V, R] = AlgEuclideJulien(P,alpha,beta)

V = 0;
R = 0;

% Matrice qui contient la valeur des fi successifs en alpha et beta 
% afin de calculer le nombre de variations pas la suite
M = zeros(2,length(P)+1); 

% Matrice qui contient les coefficients des fi successifs
F = zeros(length(P)+1,length(P)); 

F(1,:) = P;
k = numel(P)-1:-1:1;
F(2,:) = [0 P(1:end-1).*k];

M(:,1) = signeVec(polyval(F(1,:),[alpha beta]'));
M(:,2) = signeVec(polyval(F(2,:),[alpha beta]'));

i = 2;
while vectorNotNul(F(i,:))
    i = i+1;
    F(i,:) = divisionEuclidienePolynome(F(i-2,:),F(i-1,:));
    M(:,i) = signeVec(polyval(F(i,:),[alpha beta]'));
end

R =  F(i-1,:);
V = countVariations(M(1,1:i-1))-countVariations(M(2,1:i-1));
Coefficients = F(1:i-1,:);
TableauSigne = M(:,1:i-1);

end

% Fonction qui effectue la division euclidienne entre deux polynômes
function [R] = divisionEuclidienePolynome(N,D)
%% Paramètres
degNum = length(N)-1; % degré du numérateur

% le premier coefficient non nul du dénominateur ainsi que son degre
[c,degDen] = firstNonZeroCoefficient(D);


%% Initialisation
Q = 0*N; % vecteur contenant le quotient
R = N; % vecteur contenant le reste
degReste = degNum; % le degré du reste 

% Intialisation d'un vecteur de la même dimension que le numérateur
Den = 0*N;
Den(length(N)-length(D)+1:end) = D;

%% Division euclidienne : récursion
while degReste >= degDen
    facteur = R(1+degNum-degReste)/c; % facteur = premier coefficient de R non nul divisé par c
    degQ = (degReste-degDen);
    Q(end-degQ) = Q(end-degQ)+facteur;
    R = R - shiftLeft(Den*facteur, degQ);
    degReste = degReste-1;
end

R = (-1).*R; % car dans la formule de l'algorithme de Euclide c'est un - et non un +
end

% Deplace les élements d'un vecteur de nombrePlaces vers la gauche
% ex : [0 0 1 2 3] devient si nombrePlaces=1 le vecteur [0 1 2 3 0]
function v = shiftLeft(vec, nombrePlaces)
v = 0*vec;
v(1:end-nombrePlaces) = vec(1+nombrePlaces:end);
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

% Renvoie un vecteur contenant le nombre de variations de signe d'un
% vecteur
function [var] = countVariations(vec)
var = 0;
indexLastNotNul = 0;
for i = 1:length(vec)
    if vec(i)~=0
        indexLastNotNul = i;
        break;
    end
end

for i = indexLastNotNul:length(vec)
    if vec(indexLastNotNul)*vec(i)<0
        var = var +1;
        indexLastNotNul = i;
    end
end
end

% Renvoie 1 si le vecteur est nul (en prenant en compte la tolérance
% machine)
function [var] = vectorNotNul(vec)

var = 0;
tolerance = exp(-6);
for i = 1:length(vec)
    if abs(vec(i)) > tolerance
        var = 1;
        break
    end
end
end


% Renvoie un vecteur contenant les signes
function [v] = signeVec(vec)

v = vec;
tolerance = exp(-6);
for i = 1:length(vec)
    if abs(vec(i)) > tolerance
        v(i) = vec(i)/(abs(vec(i)));
    else
        v(i) = 0;
    end
end
end




