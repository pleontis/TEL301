%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all;
close all;
%A.1
%Data
T=10^(-3);
over = 10;
A=4;
Ts=T/over;
rolloff = [0 0.5 1]; 

figure(1);
hold on;
%A1
%SRRC pulse plot
for i=[1 2 3] %Create and plot srrc pulse for each roll-off factor
    [ph(i,:),t] = srrc_pulse(T,over,A,rolloff(i));
    p(i) =plot(t,ph(i,:));
end
title( 'SRRC pulses for diffrent roll-off factors' );
xlabel( 'time - t' );
ylabel( 'PH(t)' );
legend([p(1), p(2) ,p(3)], 'a = 0', 'a=0.5', 'a=1');


%A.2
%Data
Fs = 1/Ts;
Nf = 1024;
F=(-Fs/2:Fs/Nf:Fs/2-Fs/Nf);

figure(2);
hold on;

%A.2.a
%Fourier Trasnfrom
for i=[1 2 3] %Plot psd for each roll-off factor  
    FT(i,:) = fftshift(fft(ph(i,:),Nf)*Ts);
    p(i) = plot(F,abs(FT(i,:)).^2);
end

title('Plot for Power Spectral Density');
xlabel('F');
ylabel('|PH(F)|^2');
legend([p(1), p(2) ,p(3)], 'a = 0', 'a=0.5', 'a=1');

%A.2.b
figure(3);

for i = [1 2 3]     %Semilogy psd for each roll-off factor
    p(i) = semilogy(F,abs(FT(i,:)).^2);
    hold on;
end

title('Smilogy for Power Spectral Density');
xlabel('F');
ylabel('|PH(F)|^2');
legend([p(1), p(2) ,p(3)], 'a = 0', 'a=0.5', 'a=1');

%A.3.
BW1 = (1+rolloff(1))/(2*T);
BW2 = (1+rolloff(2))/(2*T);
BW3 = (1+rolloff(3))/(2*T);
fprintf('The BW for rollo-off factor a=0 is:    %d \n',BW1);
fprintf('The BW for rollo-off factor a=0.5 is:  %d \n',BW2);
fprintf('The BW for rollo-off factor a=1 is:    %d \n',BW3);

%Data
 c1=(T/10^3)*ones(1,Nf);
 c2=(T/10^5)*ones(1,Nf);
  
 figure(4);

 %Psd with semilogy and a c1
 for i = [1 2 3]    %Semilogy psd for each roll-off factor
       p(i) = semilogy(F,abs(FT(i,:)).^2);
       hold on;
 end
 semilogy(F,c1);     
 
title('Semilogy for Power Spectral Density with c1');
xlabel('F');
ylabel('|G(F)|^2');
legend([p(1), p(2) ,p(3)], 'a = 0', 'a=0.5', 'a=1');

  figure(5);
  %Psd with semilogy and a c2
 for i = [1 2 3]   
       p(i) = semilogy(F,abs(FT(i,:)).^2);
       hold on;
 end
        semilogy(F,c2);
title('Semilogy for Power Spectral Density with c2');
xlabel('F');
ylabel('|PH(F)|^2');
legend([p(1), p(2) ,p(3)], 'a = 0', 'a=0.5', 'a=1');
  
