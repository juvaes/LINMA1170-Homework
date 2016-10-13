clc; clear all; close all;

k = 0:50;
exps = 10.^(-k);

Vall = 0*k;
smallRoot = 0*k;

%% Experience pour epsilon de plus en plus petit

for i =1:length(exps)
    e = exps(i)
    coefs = [1, 1-e, e , 2*e-2*e*e,-2*e*e];
    r = roots(coefs)
    [V, R] = AlgEuclide(coefs, -10, 10);
    Vall(i)=V;
    smallRoot(i) = r(4);
end


%% Illustration
subplot(1,2,1);
plot(log(exps),Vall,'b.','MarkerSize',25); 
xlabel('Epsilon');
ylabel('Nombre de racines');
title('Graphe entre log(epsilon) et le nombre de racines réelles disctinctes trouvées')

subplot(1,2,2);
plot(log(exps),log(smallRoot),'r.','MarkerSize',25); 
xlabel('Epsilon');
ylabel('La plus petite racine trouvée');
title('Graphe log-log entre epsilon et plus petite racine trouvée')


