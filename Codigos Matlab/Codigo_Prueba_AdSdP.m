clc
close all
clear all

%% codigo de pruebas para otros codigos

CI = xlsread('Datos.xlsx','Condiciones Iniciales');
Voltajes = CI(:,3);
Angulos = CI(:,4);
Variables =  [Voltajes;Angulos];
EliminacionV = find(isnan(Voltajes(:,1)));
EliminacionO = find(isnan(Angulos(:,1)));

c = length(EliminacionV)
b = length(EliminacionO)

poti = zeros(c,1)
pot = zeros(b,1)


% for m=1:1:length(CI(:,2))
%     
%     if CI(m,2)
%     
%     else 
%         
%     if
% 
%     
% end