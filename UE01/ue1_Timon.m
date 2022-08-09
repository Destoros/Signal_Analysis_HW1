%% a&b)
m = 10;
c = 450;
d = 10;
b = 1;
a = [m d c];

H1 = tf(b,a);
for i = 1:1:2000
    w2(i) = i/10;
end
[A,B,C,D]=tf2ss(b,a);
bode(H1,w2)
grid minor;

% w = 12 bei -60 db (^= -6db)

wc = 12; % at omega about 12, you get -3db
td_lower = pi/(20*wc);  %0.0262s
td_higher = pi/(10*wc); %0.0131s, therefore i decided for Td = 0.02s


%% c)
clear all
Td = 0.02;
id = 3;
t = 0:0.001:10;
wmax = 60;
sigma = 0.001;
for i = 1:1:wmax
    F_test(i,:) = sqrt(1.9999)*sin(2*pi*i*t);                        %Creates multiple F Signals with different Frequencies
    x_out(i,:) = mass_spring_system(F_test(i,:),Td,id);         %calculates x_out with the corresponding F
    x_out_2(i,:) = x_out(i,:) + sigma*randn(1,numel(t));        %adds the noise with sigma = 0.001m 
    [ryy(i,:),my(i,:)]  = xcorr(x_out_2(i,:),'unbiased');       %acf of the noisy signal
    xmax(i) = abs(sqrt(2*(max(ryy(i,9000:11000))-sigma^2)));    %gets the amplitude of the signal without noise by 
                                                                %subtracting sigma^2 from the acf amplitude
                                                                %and then the square root of it divided by two
    h_abs(i) = 10*log(xmax(i)/max(F_test(i,:)));                %calculates the absolute value of the frequency response
    w(i) = i;                                                %calculates the corresponding omega value
    
   
    a(i,:) = fft(F_test(i,:));
    b(i,:) = fft(x_out(i,:));
    [~, inda(i)] = max(abs(a(i,:)));
    [~, indb(i)] = max(abs(b(i,:)));
    PhDiff(i) = angle(b(i,indb(i))) - angle(a(i,inda(i)));
    PhDiff(i) = PhDiff(i)*180/pi;
end
%% c)
figure(1)
subplot(2,1,1)
plot(w,h_abs)
ylabel('|H(e^{j\omega})| [dB]');
xlabel('\omega [rad/s]');
grid minor

subplot(2,1,2)
plot(w,-PhDiff)
ylabel('arg(H(e^{j\omega})[°]');
xlabel('\omega [rad/s]');
grid minor
%% d)
Td = 0.02;
id = 3;
sigma = 0.001;

F_d = 0.8*randn(2000,1);
x_d = mass_spring_system(F_d,Td,id);
x_d2 = x_d + sigma*randn(numel(F_d),1);

NDFT = 2^10;
h = hann(NDFT);
[Pff,W] = pwelch(F_d,h,NDFT/2,NDFT);
[Py1y1,W] = pwelch(x_d,h,NDFT/2,NDFT);
[Py2y2,W] = pwelch(x_d2,h,NDFT/2,NDFT);
[Pxy1,W] = cpsd(x_d,F_d,h,NDFT/2,NDFT);
[Py1x,W] = cpsd(F_d,x_d,h,NDFT/2,NDFT);
[Pxy2,W] = cpsd(x_d2,F_d,h,NDFT/2,NDFT);
[Py2x,W] = cpsd(F_d,x_d2,h,NDFT/2,NDFT);

Ha1 = sqrt(Py1y1./Pff);
Ha2 = sqrt(Py2y2./Pff);
H1xy = Pxy1./Pff;
H1yx = Py1y1./Py1x;
H2xy = Pxy2./Pff;
H2yx = Py2y2./Py2x;

figure
    plot(20*log10(abs(Ha1)))
    hold on
    plot(20*log10(abs(Ha2)))
    hold off

figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
subplot(3,1,1)
plot(W,abs(Ha1),'g','LineWidth',2)
hold on
plot(W,abs(Ha2),'r','LineWidth',2,'Linestyle','--')
legend('sqrt(Pyy/Pxx) without Noise','sqrt(Pyy/Pxx) with Noise')
axis tight
xlabel('\omega')
ylabel('|H(e^{j\omega})|')
grid on
subplot(3,1,2)
plot(W,abs(H1xy),'g','LineWidth',2)
hold on
plot(W,abs(H2xy),'r','LineWidth',2,'Linestyle','--')
legend('Pxy/Pxx without Noise','Pxy/Pxx with Noise')
axis tight
xlabel('\omega')
ylabel('|H(e^{j\omega})|')
grid on
subplot(3,1,3)
plot(W,abs(H1yx),'g','LineWidth',2)
hold on
plot(W,abs(H2yx),'r','LineWidth',2,'Linestyle','--')
legend('Pyy/Pyx without Noise','Pyy/Pyx with Noise')
axis tight
xlabel('\omega')
ylabel('|H(e^{j\omega})|')
grid on

%% e)
NDFT = 2^10;
h = hann(NDFT);
t = 0:0.001:10;
F_e = 1000*sin(2*pi*20*t);   
x_e = mass_spring_system(F_e,Td,id);
[Pxx,W] = pwelch(F_e,h,NDFT/2,NDFT);
[Pxy,W] = cpsd(x_e,F_e,h,NDFT/2,NDFT);
[Pyy,W] = pwelch(x_e,h,NDFT/2,NDFT);
coherence = (abs(Pxy).^2)./(Pxx.*Pyy);
plot(coherence)

%%
x1=sin(2*pi*2*t);
x2=sin(2*pi*2*t - pi/2);
[ccf,mccf] = xcorr(x1,x2,'biased');
[~, index_max] = max(ccf);
index_max = index_max-round(numel(ccf)/2);