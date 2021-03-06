close all; clc; clear all;

LASTN = maxNumCompThreads(1);

%% Donn�es

alpha = -1000;
beta = 1000;

% vecteur contenant le nombre de racines que l'on va introduire
numberRoots = 500:500:4000;
timeExecutionCas1 = 0*numberRoots;
timeExecutionCas2 = 0*numberRoots;

%% Exp�rience : 2 cas diff�rents

for i = 1:length(numberRoots)
    n = numberRoots(i)
    % Cas 1 : on prend un des coefficients entiers au hasar entre 0 et 10
    coefficients = randi(10,1,numberRoots(i)); 
    tic;
    [V,R] = AlgEuclide(coefficients,alpha,beta);
    timeExecutionCas1(i) = toc;
    
    % Cas 2 : on prend des coefficients �gaux � 1 (cas m�chant)
    coefficients = ones(1,numberRoots(i)); 
    tic;
    [V,R] = AlgEuclide(coefficients,alpha,beta);
    timeExecutionCas2(i) = toc;
    
end
 
%% Calcul de la r�gression lin�naire
regressionLinCas1 = polyfit(log(numberRoots),log(timeExecutionCas1),1)
regressionLinCas2 = polyfit(log(numberRoots),log(timeExecutionCas2),1)

x = linspace(log(numberRoots(1)),log(numberRoots(end)), 10);
yCas1 = polyval(regressionLinCas1,x);
yCas2 = polyval(regressionLinCas2,x);

%% Illustration de la r�gression lin�aire
subplot(1,2,1);
plot(log(numberRoots),log(timeExecutionCas1),'r.','MarkerSize',35); hold on;
xlabel('Nombre de racines');
ylabel('Temps d excution (s)');
title('Graphe log-log entre le nombre de racines et le temps d ex�cution')
plot(x,yCas1);

subplot(1,2,2);
plot(log(numberRoots),log(timeExecutionCas2),'r.','MarkerSize',35); hold on;
xlabel('Nombre de racines');
ylabel('Temps d excution (s)');
title('Graphe log-log entre le nombre de racines et le temps d ex�cution')
plot(x,yCas2);