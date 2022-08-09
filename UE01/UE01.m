close all
clear all
clc

%if the folder Figures does not exists, create it
mkdir('Figures')    

    



% Reindl Hannes
% 11.11.2019
% 01532129
% UE01


%my ID number
id = 25;

% (A) ---------------------------------------------------------------------
% derive state space representation
m = 10;
c = 450;
d = 10;


H = tf([1/m],[1 d/m c/m]);


%check if correct with inbuild func
% [b,a] = tfdata(H,'v');
% [A_test, b_test, c_T_test, d_test] = tf2ss(b,a);
%OR WITH
% sys_2 = ss(H); 

%State Space representation derived through "1. Normalform" and Transfer
%Function
A = [  0    1;
     -c/m -d/m];
 
b = [0; 1];
c_T = [1/m 0];
d_through = 0;

%convert to build in Matlab State Space representation
sys = ss(A,b,c_T,d_through);

%Plot the frequency response
figure
    bodeplot(sys)
    grid on
    
    saveas(gcf,'Figures/analytical_freqeuncy_response', 'epsc')
    
%looks like a lowpass filter with a increase in the magnitude around the
%resonance frequency of the system
%its a second order 'lowpass', because the decrease in magnitude is
%-40dB/decade and the phase shift goes to -180°
%The pass region of the frequency response is already decreased by about -50dB


% (B) ---------------------------------------------------------------------
% determine bandwidth of system

% OmegaRes = imag((-d + sqrt(d^2-4*c*m))/(2*m)); %formula from paper
OmegaG = 10; %read from plot;

Td = pi/(20*OmegaG); %given heuristic






% (C) ---------------------------------------------------------------------
% estimate frequency response through steady state



t = 0:Td:100; %Td from (B); 100sec to have enough samples
F_hat = sqrt(1.9); %Expectation of F^2 should be less than 1



ii = 1;
f_max = 1/(2*Td) ;
f_plot = logspace(-1,log10(f_max/4),200);

H_abs = zeros(1,length(f_plot));
H_phi = zeros(1,length(f_plot));

%calculate input to output realtion for several different frequencies
for f = f_plot
    
 
    F = F_hat * sin(2*pi*f*t);% Force Vector   


    y=mass_spring_system(F,Td,id);  
    
%     figure
%         plot(t,y)
        %seen from plot: after 20 seconds, the steady state solution
        %settles 
        
    index_steady_state = ceil(20/Td);
    
    %create noise to corrupt y
    noise = randn(length(t),1)/1000; %variance is equal to AC Power
    sigma_squared = var(noise);
    
    %add measurement noise
    y = y+noise; 
    
    y_steady_state = y(index_steady_state:end);
    
    %compute acf to calculate the peak value
    [ryy, myy] = xcorr(y_steady_state,'biased'); %acf has some immunity to noise,
    %since we only look at ryy[0] the biased doesnt change the value
    
%     figure
%         plot(myy,ryy)
    
    S_plus_N = ryy(myy==0); %ryy[0] is equal to power of signal + power of noise
      
    if S_plus_N - sigma_squared < 0 %can only occur due to low signal power and numerical errors
        %if this occurs ,keep the same value as in the previous loop
        disp('here')
    else
        S = S_plus_N - sigma_squared;  %substract noise power(variance is equal to AC Power)
    end
    
    H_abs(ii) = sqrt(2)*sqrt(S)/F_hat; %sqrt(S) = U_eff; sqrt(2) to get peak value

        
    F_steady_state = F(index_steady_state:end);
    [ryF,myF] = xcorr(y_steady_state,F_steady_state,'biased');
    %get phase shift; first peak is the shortest phase shift


    
    [~, index_max] = max(ryF); %using the biased estimater the closest peak should be the maximum value
  

    T = 1/f;
    H_phi(ii)  = -sign(index_max)*mod(myF(index_max)*Td,T)*360/T;

    
    
    ii = ii + 1;
%     
%     figure
%         plot(t,50*y)
%         hold on
%         plot(t,F)
%         hold off
%         title('Input signal F(t) and output signal y(t)')
%         xlabel('t')
%         legend('100 \cdot y(t)', 'F(t)')
%         xlim([0 5])
%         grid on
%         
%     figure
%         plot(myF,ryF)
%         title('KKF of y(t) and F(t)')
%         xlabel('myF')
%         legend('ryF[m]')
%         grid on

  

end


 H_dB = 20*log10(H_abs);
        
figure
    subplot(2,1,1)
    semilogx(f_plot,H_dB)
    title('Estimated frequency response')
    ylim([-100 -20])
    grid on
    ylabel('Magnitude in dB')
    
    subplot(2,1,2)
    semilogx(f_plot,H_phi)
    grid on
    ylabel('Phase in deg')
    
    xlabel('Frequecy')

    saveas(gcf,'Figures/steady_state_freq_resp', 'epsc')
%Mithilfe der bekannten steady state Lösung kann der Betrag und die Phase bestimmt werden.
%Der Betrag wird aus den Amplituden unterschied vor und nach dem
%mass_spring_system() berechnet. Die Ampltiduen werden durch die AKF berechnet
%Die Phase wird mithilfe der KKF berechnet. Die
%steady state Lösung schwingt genau mit der gleichen Frequenz wie die Kraft
%F, dadurch ist nur ein Phasenunterschied zwischen Ein- und Ausgang,
%welcher einen Shift der beiden Signale entspricht. Mithilfe der KKF(biased) wird
%dann der größte Peak gesucht, welcher aufgrund der biased Methode auch der
%nächstegelegne von shift = 0 ist. Daduch erhält man die Anzahl der
%geshifteten Samples und kann mithilfe der Sampleperiode Td auf die Zeit
%Zurückrechen. Von dieser Zeit und der Periodendauer kann dann die
%Phasenverchiebung berechnet werden.
%Der modulo operator dient nur zur Absicherung, falls mal der größte Peak (z.B.
%aufgrund von Rauschen) nicht der nächstgelegen ist. Der modulo operator
%"subtrahiert" solange die Periodendauer des Signal T, bis es innnerhalb
%einer Periode ist.

%plot beschreibung:
%Der aus der Transferfunktion mit dem vorgegeben Werten
%und der geschätzte Frequenzgang ähnlich sich auf den ersten Blick ziemlich
%gut. Die Grenz- bzw Resonanzfrequenz sind natürlich unterschiedlich, da
%die gegebene Werte für c,d und m nicht mit denen des mass_spring_system()
%mit id = 25 übereinstimmen. Die -90° Phasendrehung bei der Resonanzfrequenz 
%passt auch gut.
%Auch fällt auf, das im geschätzen Frequnezgang die Phasenverschiebung 
% unter -180° geht und mit fortlaufenden f immer mehr negativ wird. Wie
% dies Zustanden kommt kann ich mir nicht erklären.
%-----------------------------------HIER NOCHMAL DENKEN HANNES------------


% (D) ---------------------------------------------------------------------
% estimate frequency response thorugh white noise

for ii = 1:2
    x = randn(length(t),1); %variance for randn is 1


    y=mass_spring_system(x,Td,id);

    %first time dont corrupt the output signal with noise
    if ii == 1
       y = y; %no measurement noise
       text_title = ['estimated frequency response w/o measurement noise'];
    else
          y = y +  randn(length(t),1)/1000; %add measurement noise with sigma = 0.001; %create new noise, so that the input noise and measurement noise are not correlated
%         y = y +  x/1000; %test with correlated noise
        text_title = ['estimated frequency response with measurement noise'];
    end

%     % Correlation using unbiased estimator
%     [rxx,mx]  = xcorr(x,'unbiased'); %why unbiased and not biased?
%     [ryy,my]  = xcorr(y,'unbiased');
%     [ryx,myx] = xcorr(y,x,'unbiased');
% 
%     
%     scale_fac = 1000;
% 
%     figure, hold on, set(gca,'FontSize',26),set(gcf,'Color','White');
%         title('Auto- and Crosscorrelation')
%         plot(mx,rxx,'b','LineWidth',2)
%         plot(my,ryy*scale_fac,'r','LineWidth',2)
%         plot(myx,ryx*scale_fac,'g','LineWidth',2)
%         axis tight
%         xlabel('m')
%         legend('r_{xx}','r_{yy} scaled','r_{yx} scaled')
%         grid on


    NDFT = 2^10; %tried several NDFT values until the plot looked quite well
    h = hann(NDFT);

    [Pxx,W] = pwelch(x,h,NDFT/2,NDFT); %using pwelch like in the example Matlab files
    [Pyy,W] = pwelch(y,h,NDFT/2,NDFT);
    [Pxy,W] = cpsd(y,x,h,NDFT/2,NDFT);


%     H_est_2 = abs(Pxy./Pyy); %absolut value of H
%     H_est_2_dB = 20*log10(H_est_2);
    
    H_est_2 = sqrt(Pyy./Pxx); %super formula aka Wiener Lee eq
    H_est_dB = 20*log10(H_est_2);

    figure
        semilogx(W, H_est_dB)
        grid on
        title(text_title)
        xlabel('normalized frequency to \pi')
        ylabel('|H(e^{j \omega})|')
        ylim([-100 -20])
        
        saveas(gcf,['Figures/' strrep(strrep(text_title,' ', '_'),'w/o', 'without')], 'epsc')
        
% (E) ---------------------------------------------------------------------
% coherence function

gamma(:,ii) = abs(Pxy).^2./(Pxx .* Pyy);
        
end



gamma_mean_1 = mean(gamma(:,1))
gamma_mean_2 = mean(gamma(:,2))

figure
%     plot(W,gamma(:,1),'LineWidth', 1.5)
    plot(W,gamma(:,1),'LineWidth', 1.5)
    hold on
%     plot(W,gamma(:,2),'LineWidth', 1.5)
    plot(W,gamma(:,2),'LineWidth', 1.5)
    title('Coherence')
    grid on
    hold off
    legend('gamma w/o measurement noise', 'gamma with measurement noise')
    xlabel('normalized frequency to \pi')
    ylabel('coherence \gamma')
    
    saveas(gcf,'Figures/coherence', 'epsc')
 
    
disp('done')
%0 <= gamma <= 1
%für lineares system: gamma = 1
%im ersten durchlauf, also ohne measurement noise ist gamma ziemlich
%nah an 1 dran über das gesammte Spektrum.(stimmt nicht für Randwerte, aber
%ich geh mal davon aus, das hat mit dem begrenzten Signal und Windowing zu
%tun)

%im zweiten durchlauf variiert das gamma zwischen 1 und 0, wobei es eher im
%unteren Bereich ist



F_test = 0.5*t;
y=mass_spring_system(F_test,Td,id);



[Pxx,W] = pwelch(F_test,h,NDFT/2,NDFT); %using pwelch like in the example Matlab files
[Pyy,W] = pwelch(y,h,NDFT/2,NDFT);
[Pxy,W] = cpsd(y,F_test,h,NDFT/2,NDFT);

gamma_test = abs(Pxy).^2./(Pxx .* Pyy);

figure
    plot(W,gamma_test)
    title('Coherence')
    grid on
    hold off
    legend('Coherence for linear F')
    xlabel('normalized frequency to \pi')
    ylabel('coherence \gamma')
    
    saveas(gcf,'Figures/coherence_test', 'epsc')
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

%w = 2*pi*f; f= w/2/pi; fs > 2 * f where f is the highest occuring
%frequency; Td = 1/fs; 
%fs > 2 * w/(2*pi)
%Td < 2*pi/(2*w)
%Td < pi/w; 





% [z,p,k] = zpkdata(H, 'v'); %the v causes the given data in z,p,k to be normal arrays and not cell arrays
% OmegaG can be retrieved from the poles

% figure
% [b,a] = tfdata(H, 'v')
%  freqz(b,a,512);

% figure 
% pzmap(H)




%derive OmegaG from plot

% w = logspace(0,1,10000);
% [mag,phase,wout] = bode(sys,w); %plots frequency response
% OmegaG = w(max(mag) == mag) %Not correct