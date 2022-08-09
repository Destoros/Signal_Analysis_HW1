function [ y ] = func_Quantizer( Q, x )
%
% Quantizer
%
%  Q: Quantizer
%  x: Input signal
%
% Neumayer 2013


y = zeros(size(x));

%Clipping
idx1 = find(x>=Q.u_in(end)); 
x(idx1) = Q.delta;

idx2 = find(x<=Q.u_in(1)); 
x(idx2) = Q.delta;


for ii = 1: max(size(x))
     
   idx = find(Q.u_in<x(ii));
   y(ii) = Q.u_out(idx(end)+1);
 
    
end
y(idx1) = Q.u_out(end);
y(idx2) = Q.u_out(1);

end

