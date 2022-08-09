% SA Problem Class
%
% Demo Example for the function func_ADC
%
%
% Neumayer 2017
clear all, close all, clc


% The following lines show the implementation of the function func_ADC

% Nr Samples
N = 1000;


% Create ADC Object
X_m = 2; B   = 4;
ADC = func_createQuantizer(X_m,B); 

% Plot quantizer curve
func_plotQuantizer( ADC ,1 ,'b')


% Sampling Frequency
ADC.f_S = 1E3;

%T Jitter
ADC.Tjitter = 1E-6;   % T Jitter varies with your ID number


% Sinusoidal Signal
sig = @(t,PAR) PAR.A*sin(2*pi*PAR.f*t+PAR.phi);

PAR.A = 1;      %Amplitude
PAR.f = 39;     % Frequency 
PAR.phi = 0;    % Phase


t = (0:N-1)/ADC.f_S;
x = sig(t,PAR);

tsuh = t + ADC.Tjitter * (2*rand(size(t))-1)  ;

xSH = sig(tsuh,PAR);

[ xQ   ] = func_Quantizer( ADC, x );
[ xSHQ ] = func_Quantizer( ADC, xSH );

% How to run func_ADC
id = 52; % your ID
[x,xSH,xQ,xSHQ] = func_ADC(sig,PAR,N,id);


help func_ADC

