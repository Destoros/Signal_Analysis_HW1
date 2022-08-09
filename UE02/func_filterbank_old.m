function [x1,t1,x2,t2,x3,t3,x4,t4] = func_filterbank_old(x,f_S,h0,h1)

%fitlerbank structure
x1 = upfirdn(x,h1,1,2); %x1 = x1(1:N/2); %delete entries which assumed x[n] = 0 for n > length(x)
x1_5 = upfirdn(x,h0,1,2);% x1_5 = x1_5(1:N/2);

x2 = upfirdn(x1_5,h1,1,2);% x2 = x2(1:N/4);
x2_5 = upfirdn(x1_5,h0,1,2);% x2_5 = x2_5(1:N/4);

x3 = upfirdn(x2_5,h1,1,2);% x3 = x3(1:N/8);
x4 = upfirdn(x2_5,h0,1,2);% x4 = x4(1:N/8);


%create new time vectors
T_S = 1/f_S;

t1 = T_S*2 * (1:length(x1)); %1 times downsampled with L = 2; 2^1 = 2
t2 = T_S*4 * (1:length(x2)); %2 times downsampled with L = 2; 2^2 = 4
t3 = T_S*8 * (1:length(x3)); %3 times downsampled with L = 2; 2^3 = 8
t4 = T_S*8 * (1:length(x4)); %3 times downsampled with L = 2; 2^3 = 8

end