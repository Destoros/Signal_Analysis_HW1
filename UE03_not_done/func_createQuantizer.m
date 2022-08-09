function Q = func_createQuantizer(X_m,B)
%
% func_creatQuantizer creates an ideal bipolar quantizer 
% 
% X_m ... Full Scale Range
%  B  ... Number of Bits
%
% Neumayer 28.7.2013


Q.X_m = X_m;
Q.B   = B;

N = 2^B;
delta = X_m/N;

Q.delta = delta;

%Quantization
Q.u_in  = [-X_m/2 + delta/2: delta: X_m/2- 1.5*delta];
Q.u_out = [-X_m/2: delta: X_m/2 - delta]; 

%Q.u_in  = [-X_m/2 - delta/2: delta: X_m/2 - 1.5*delta/2];
%Q.u_out = [-X_m/2: delta: X_m/2 - delta]; 