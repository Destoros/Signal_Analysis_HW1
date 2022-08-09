% SA Problem Class
%
% Demo Example on the use of mass_spring_system
%
%
% Neumayer 2019
clear all, close all, clc

% Force Vector
F = [ones(1,1);zeros(100,1)];


Td = 0.1;

id = 1;

y=mass_spring_system(F,Td,id);

t = [0:1:length(y)-1]*Td;


figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
plot(t,y,'r','LineWidth',2), grid on
xlabel('t (s)')
ylabel('y(t) (m)')
