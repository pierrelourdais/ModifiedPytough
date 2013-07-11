close all
clear all
path = '/Users/Pierre/Desktop/CO2';
load ([path, '/co2data3.mat']);

Reynolds_limit=10.0;
porosity=.41;
A=pi*(.045^2);
d_50=0.000105;
%Dynamic viscosity is given by co2data3.mat <--- M

Qmax= (Reynolds_limit*porosity*A*M)/(d_50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Range plot of T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tplot=zeros(20,15);
for i=1:15
Tplot(:,i)=5*i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Range plot of P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pplot=zeros(20,15);
for j=1:20
    Pplot(j,:)=5*j;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Difference calcul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qmaxplot=zeros(20,15);

for k=1:15
    for l=1:20
       Qmaxplot(l,k)= interp2(T,P,Qmax,Tplot(l,k),Pplot(l,k)); 
    end
end



figure
[c,h] = contour(Tplot,Pplot,Qmaxplot);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
title('Qmax to get Darcy law (kg/s');
ylabel('Pressure (bar)');
xlabel('Temperature (C)');