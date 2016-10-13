close all; clc;

% On impose à Matlab d'utiliser un seul thread
LASTN = maxNumCompThreads(1);

%% Données

% vecteur contenant le nombre de racines que l'on va introduire
numberRoots = 500:500:4000;
timeExecutionCas1 = 0*numberRoots;
timeExecutionCas2 = 0*numberRoots;

%% Expérience : 2 cas différents

for i = 1:length(numberRoots)
    
    % Cas 1 : on prend un des coefficients entiers au hasar entre 0 et 10
    coefficients = randi(10,1,numberRoots(i)); 
    tic;
    roots(coefficients);
    timeExecutionCas1(i) = toc;
    
    % Cas 2 : on prend des coefficients égaux à 1 (cas méchant)
    coefficients = ones(1,numberRoots(i)); 
    tic;
    roots(coefficients);
    timeExecutionCas2(i) = toc;
    
end
 
%% Calcul de la régression linénaire
regressionLinCas1 = polyfit(log(numberRoots),log(timeExecutionCas1),1)
regressionLinCas2 = polyfit(log(numberRoots),log(timeExecutionCas2),1)

x = linspace(log(numberRoots(1)),log(numberRoots(end)), 10);
yCas1 = polyval(regressionLinCas1,x);
yCas2 = polyval(regressionLinCas2,x);

%% Illustration de la régression linéaire
subplot(1,2,1);
plot(log(numberRoots),log(timeExecutionCas1),'r.','MarkerSize',35); hold on;
xlabel('Nombre de racines');
ylabel('Temps d excution (s)');
title('Graphe log-log entre le nombre de racines et le temps d exécution')
plot(x,yCas1);

subplot(1,2,2);
plot(log(numberRoots),log(timeExecutionCas2),'r.','MarkerSize',35); hold on;
xlabel('Nombre de racines');
ylabel('Temps d excution (s)');
title('Graphe log-log entre le nombre de racines et le temps d exécution')
plot(x,yCas2);

