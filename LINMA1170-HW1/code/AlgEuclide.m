function [V, R] = AlgEuclide(P,alpha,beta)

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
    R = makeTrueZero(divisionEuclidienePolynome(F(i-2,:),F(i-1,:)));
    F(i,end-length(R)+1:end) = R;
    M(:,i) = signeVec(polyval(F(i,:),[alpha beta]'));
end

R =  F(i-1,:);
V = countVariations(M(1,1:i-1))-countVariations(M(2,1:i-1));


% Va = countVariations(M(1,1:i-1));
% Vb = countVariations(M(2,1:i-1));
% Coefficients = F(1:i-1,:)
% TableauSigne = M(:,1:i-1);

end

% Fonction qui effectue la division euclidienne entre deux polynômes
function [R] = divisionEuclidienePolynome(N,D)

%% Tests preliminaires
NwithOutZero=WithoutZero(N);
DwithOutZero=WithoutZero(D);
if DwithOutZero(1)==0 || NwithOutZero(1)==0 %Les vecteurs sont vides
    error('error: Vecteurs vides ou division par 0');
    R=0;
elseif length(NwithOutZero)<length(DwithOutZero) %Si le polynôme N a un degré plus faible que le polynôme D, alors la division euclidienne est aussi impossible.
    error('error: le polynome N a un degré plus faible que le polynome D : le quotient est nul');
    R=NwithOutZero;
else
    k=length(NwithOutZero)-length(DwithOutZero); 
    newN=NwithOutZero;%Création d'un nouveau polynôme N dans le but de ne pas modifier le polynôme N déjà existant
    % La division euclidienne.
    for i = 1:k+1
        coeff = newN(i)/DwithOutZero(1); %Calcul du coefficient
        newD = zeros(length(NwithOutZero),1); %Création d'un nouveau polynôme D dans le but de ne pas modifier le polynôme D déjà existant
        %Comme le degré du polynôme D varie à chaque étape de la division euclidienne, nous allons
        % en créer 1 nouveau à chaque étape
        newD(i:length(newD)+i-k-1) = DwithOutZero(:);%On va changer les exposants des termes de D (principe de la division euclidienne).
        newN = newN-(newD.* coeff)';
    end
    R=WithoutZero(newN);%Suppression des termes inutiles.
end

% Pour eviter de propager les erreurs numeriques
if length(R)>length(DwithOutZero)
    R = R(end-length(DwithOutZero)+1:end);
end
R = (-1)*makeTrueZero(R);
end

function newVec = makeTrueZero(vec)
newVec = vec;
for j = 1:length(vec)
    if abs(vec(j))<10^(-8)
        newVec(j)=0;
    end
end
end

%Fonction servant à éliminer les termes nuls qui doivent être
%supprimés dans un vecteur vec, c'est-à-dire tous les coefficients nuls devant le premier
%coefficient (non-nul) de l'indice le plus petit du vecteur. Si le vecteur vec est
%vide ou n'est composé que de 0, le vecteur retourné sera [0].
function [newVec] = WithoutZero(vec)
if isempty(vec)
    newVec=[0];
else
    i=1;
    while vec(i)==0 && i<length(vec)
        i=i+1;
    end
    for j=i:length(vec)
        newVec(j-i+1)=vec(j);
    end
end
end

% Fonction renvoyant un vecteur contenant le nombre de variations de signe
function [var] = countVariations(vec)
var = 0;
indexLastNotNul = 0;
for i = 1:length(vec)
    if vec(i)~=0
        indexLastNotNul = i;
        break;
    end
end

for i = indexLastNotNul:length(vec)-1
    if vec(indexLastNotNul)*vec(i+1)<0
        var = var+1;
    end
    
    if vec(i)~=0
        indexLastNotNul = i+1;
    end
end
end

% Fonction qui renvoie true si le vecteur est non nul
function [var] = vectorNotNul(vec)
var = 0;
for i = 1:length(vec)
    if abs(vec(i)) > 0
        var = 1;
        break
    end
end
end


% Renvoie un vecteur contenant les signes d'un vecteur [-5 0 1 2] --> [-1 0 1 1]
function [v] = signeVec(vec)
v = 0*vec;
for i = 1:length(vec)
    if vec(i) > 0
        v(i) = 1;
    elseif vec(i) < 0
        v(i) = -1;
    end
end
end




