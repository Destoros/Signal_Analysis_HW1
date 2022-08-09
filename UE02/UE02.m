close all
clear all
clc

%if the folder Figures does not exists, create it
mkdir('Figures')    


% Reindl Hannes
% 11.11.2019
% 01532129
% UE02

%ID = 25

%--------------------------------------------------------------------------
%(A)

H0 = Lowpass_filter(); %designed with fda tool; equiripple filter
H1 = Highpass_filter();

[h0_plot, w0] = freqz(H0);
[h1_plot, w1] = freqz(H1);

h0_dB = 20*log(abs(h0_plot));
h1_dB = 20*log(abs(h1_plot));

figure
    plot(w0,h0_dB)
    hold on
    plot(w1,h1_dB)
    grid on
    title('(A) Frequency response of H_0(z) and H_1(z)')
    xlabel('frequency normalized to \pi')
    ylabel('$|H|$ in dB','Interpreter','latex')
    legend([{'$|H_0(e^{j \omega})|$'}, {'$|H_1(e^{j \omega})|$'}],'Interpreter','latex')
    xlim([0, pi])
    
    saveas(gcf,'Figures/a_frequency_response', 'epsc')
    



%--------------------------------------------------------------------------
%(B)
%siehe inkscape plot

%--------------------------------------------------------------------------
%(C)
%see function  func_filterbank(x,f_S,h0,h1)




%--------------------------------------------------------------------------
% (D)

f_S = 2000;
f_max = f_S/2;

%Assuming the filters perfect filters: (which they are not)
%could be computed but not worth the effort
Bandwidth_x1 = [f_max/2, f_max]; %upper and lower frequency boundary for signal x1
Bandwidth_x2 = [f_max/4, f_max/2];
Bandwidth_x3 = [f_max/8, f_max/4];
Bandwidth_x4 = [0, f_max/8];


Td = 1/f_S;
f_max = f_S/2;
f1 = mean(Bandwidth_x1); %put the frequency of the sinusoidal in the middle of the bandwidth
f2 = mean(Bandwidth_x2);
f3 = mean(Bandwidth_x3);
f4 = mean(Bandwidth_x4);

f1 = 900; %for plot, to describe the mirroring

T_max = 1/min([f1,f2,f3,f4]);

t = 0:Td:10*T_max; %show 2 full periods of the lowest frequency


x = 1*sin(2*pi*f1*t) + 1*sin(2*pi*f2*t) + 1*sin(2*pi*f3*t) + 1*sin(2*pi*f4*t);
figure
    spectrogram(x,256,0,256,f_S)
    title('(D) Spectrogramm of sinusoidals')
    
     saveas(gcf,'Figures/d_spectromgramm', 'epsc')

h0_num = H0.Numerator;
h1_num = H1.Numerator;

[x1,t1,x2,t2,x3,t3,x4,t4] = func_filterbank(x,f_S,h0_num,h1_num);

X = {x1, x2, x3, x4};
T = {t1, t2, t3, t4};


figure
for ii = 1:length(X)
    
    subplot(length(X),1,ii)
    plot(T{ii},X{ii})
    ylabel(['x_' num2str(ii) '[n]'])
    
    if ii == 1; title('(D) filter bank outputs for different sinusoidal'); end
    
    if ii == length(X); xlabel('time'); end
    
end
    saveas(gcf,'Figures/d_filter_bank_output', 'epsc')



%--------------------------------------------------------------------------
% (E)
UE2sig = load('UE2sig');
UE2sig = UE2sig.x;
f_S = 2000;

figure
    spectrogram(UE2sig,256,0,256,f_S)
    title('(E) Spectrogramm of given signal')
    
    saveas(gcf,'Figures/e_spectromgramm', 'epsc')

[x1,t1,x2,t2,x3,t3,x4,t4] = func_filterbank(UE2sig,f_S,h0_num,h1_num);

X = {x1, x2, x3, x4};
T = {t1, t2, t3, t4};

figure
for ii = 1:length(X)
    
    subplot(length(X),1,ii)
    plot(T{ii},X{ii})
    ylabel(['x_' num2str(ii) '[n]'])
    
    if ii == 1; title('(E) filter bank outputs for given signal'); end
    
    if ii == length(X); xlabel('time'); end
    
end

    saveas(gcf,'Figures/e_filter_bank_output', 'epsc')







