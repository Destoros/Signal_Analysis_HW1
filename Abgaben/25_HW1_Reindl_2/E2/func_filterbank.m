function [x1,t1,x2,t2,x3,t3,x4,t4] = func_filterbank(x,f_S,h0,h1)


M0 = length(h0);
M1 = length(h1);

%fitlerbank structure
%the downsampling of the highpass filtered signal is not trivial. For the
%change of the frequency spectrum look at the attached PDF: Downsampling_change_frequency_domain
x1 = upfirdn(x,h1,1,2);
x1 = x1(ceil(M1/4):end-floor(M1/4)); %delete entries which assumed x[n] = 0 

x1_5 = upfirdn(x,h0,1,2);% x1_5 = x1_5(1:N/2);
x1_5 = x1_5(ceil(M0/4):(end - floor(M0/4))); %delete entries which assumed x[n] = 0 



x2 = upfirdn(x1_5,h1,1,2);% x2 = x2(1:N/4);
x2 = x2(ceil(M1/4):(end - floor(M1/4))); %delete entries which assumed x[n] = 0 

x2_5 = upfirdn(x1_5,h0,1,2);% x2_5 = x2_5(1:N/4);
x2_5 = x2_5(ceil(M0/4):(end - floor(M0/4))); %delete entries which assumed x[n] = 0 



x3 = upfirdn(x2_5,h1,1,2);% x3 = x3(1:N/8);
x3 = x3(ceil(M1/4):(end - floor(M1/4))); %delete entries which assumed x[n] = 0 

x4 = upfirdn(x2_5,h0,1,2);% x4 = x4(1:N/8);
x4 = x4(ceil(M0/4):(end - floor(M0/4))); %delete entries which assumed x[n] = 0 



%create new time vectors
T_S = 1/f_S;

t1 = T_S*2 * (1:length(x1)); %1 times downsampled with L = 2; 2^1 = 2
t2 = T_S*4 * (1:length(x2)); %2 times downsampled with L = 2; 2^2 = 4
t3 = T_S*8 * (1:length(x3)); %3 times downsampled with L = 2; 2^3 = 8
t4 = T_S*8 * (1:length(x4)); %3 times downsampled with L = 2; 2^3 = 8



end