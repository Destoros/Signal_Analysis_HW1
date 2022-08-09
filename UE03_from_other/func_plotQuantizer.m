function   func_plotQuantizer( Q ,fignr ,colstr)
% Plots the quantization curve of the ADC Object
% 
%
% fignr = figure number
% colstr = Color
%
% Neumayer 2013

u_in  = Q.u_in;
u_out = Q.u_out;
delta = Q.delta;

figure(fignr), hold on, set(gcf,'Color','White');
set(gca,'FontSize',26)

plot([u_in(1)-2*delta, u_in(1)],[u_out(1),u_out(1)],colstr,'LineWidth',2)
for ii = 1:length(u_in)
 % plot([u_in(ii), u_in(ii+1)],[u_out(ii), u_out(ii)])
  plot([u_in(ii), u_in(ii)],[u_out(ii), u_out(ii+1)],colstr,'LineWidth',2)  
end
 plot([u_in(end), u_in(end)+2*delta],[u_out(end),u_out(end)],colstr,'LineWidth',2)

 for ii = 1:length(u_in)-1
 % plot([u_in(ii), u_in(ii+1)],[u_out(ii), u_out(ii)])
  plot([u_in(ii), u_in(ii+1)],[u_out(ii+1), u_out(ii+1)],colstr,'LineWidth',2)  
end
 
 
grid on
h = axis;
h = axis *1.2;
axis(h);
axis equal
