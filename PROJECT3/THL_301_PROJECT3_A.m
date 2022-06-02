%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all;
close all;
%% Specifications
N = 100;
Nf = 2048;
T = 10^(-3);
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 4;
a = 0.5;
Fo = 2000;

%% TASK 1
b = (sign(randn(N,3))+1)/2;
%% Task 2
Xn=bits_to_PSK_8(b);
XI=Xn(:,1);
XQ=Xn(:,2);
%% Task 3
%Create SRRC pulse 
[ph,t] = srrc_pulse(T,over,A,a);
F = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;
%Create Xd signals for usage in Convolution
Xdi = Fs*upsample(XI,over);
t_Xdi =(0:Ts:N/(1/T)-Ts);
Xdq = Fs*upsample(XQ,over);
t_Xdq = (0:Ts:N/(1/T)-Ts);
%Compute Convolution and create time axis. Also find length of Conv.
t_Xti = (t(1)+t_Xdi(1):Ts:t(end)+t_Xdi(end));
t_Xtq = (t(1)+t_Xdq(1):Ts:t(end)+t_Xdq(end));
Xti = conv(ph,Xdi)*Ts;
Xtq = conv(ph,Xdq)*Ts;
Xti_total=length(t_Xti)*Ts;
Xtq_total = length(t_Xtq)*Ts;
%Fourier Transform and compute Periodogram for both signals 
Fxi=fftshift(fft(Xti,Nf)*Ts);
Pxi=(abs(Fxi).^2)/Xti_total;
Fxq=fftshift(fft(Xtq,Nf)*Ts);
Pxq=(abs(Fxq).^2)/Xtq_total;
%Show results for each signal
figure();
subplot(3,1,1);
plot(t_Xdi,Xdi);
title('Xi[n]');
subplot(3,1,2);
plot(t_Xti,Xti);
title('Xi,n after filtering with SRRC pulse');
xlabel('t');
subplot(3,1,3);
plot(F,Pxi);
title('Periodogram of Xi,n after filtering with SRRC pulse'); 
ylabel('Px(F)');
xlabel('Frequency (Hz)');

figure();
subplot(3,1,1);
plot(t_Xdq,Xdq);
title('Xq[n]');
subplot(3,1,2);
plot(t_Xtq,Xtq);
title('Xq,n after filtering with SRRC pulse');
xlabel('t');
subplot(3,1,3);
plot(F,Pxq);
title('Periodogram of Xq,n after filtering with SRRC pulse'); 
ylabel('Px(F)');
xlabel('Frequency (Hz)');
%% TASK 4
%Create XI(t) and XQ(t) signals
XTi = 2*(Xti).*(cos(2*pi*Fo*transpose(t_Xti)));
Xti_total = length(t_Xti)*Ts;

XTq = -2*(Xtq).*(sin(2*pi*Fo*transpose(t_Xtq)));
Xtq_total = length(t_Xtq)*Ts;
%Fourier Transform and compute Periodogram for both signals
Fxi = fftshift(fft(XTi,Nf)*Ts);
Pxi = (abs(Fxi).^2)/Xti_total;

Fxq = fftshift(fft(XTq,Nf)*Ts);
Pxq = (abs(Fxq).^2)/Xtq_total;
%Show results for each signal
figure();
subplot(2,1,1);
plot(t_Xti,XTi);
title('XI(t)');
xlabel('t');
subplot(2,1,2);
plot(F,Pxi);
title('Periodogram for XI(t)'); 
ylabel('Px(F)');
xlabel('Frequency (Hz)');

figure();
subplot(2,1,1);
plot(t_Xtq,XTq);
title('XQ(t)')
xlabel('t');
subplot(2,1,2);
plot(F,Pxq);
title('Periodogram for XQ(t)'); 
ylabel('Px(F)');
xlabel('Frequency (Hz)');

%% TASK 5
%Create channel input signal
Xt = XTi+XTq;
T_total = length(t_Xtq)*Ts;
%Fourier Transform and compute Periodogram
Fxt = fftshift(fft(Xt,Nf)*Ts);
Pxt = (abs(Fxt).^2)/T_total;

figure();
subplot(2,1,1);
plot(t_Xtq,Xt);
title('X(t)=XI(t)+XQ(t)');
xlabel('t');
subplot(2,1,2);
plot(F,Pxt);
title('Periodogram for X(t)');
ylabel('Px(F)');
xlabel('Frequency (Hz)');

%% TASK 7
%Compute variance and create Gaussian Noise
SNR = 20;
s2w = 1/(Ts*(10^(SNR/10)));
Wt = sqrt(s2w).*randn(length(Xt),1);
%Create Y(t) signal
Yt = Xt+Wt;
%Show results
figure();
plot(t_Xtq,Xt,'red');
hold on;
plot(t_Xtq,Yt,'blue');
hold off;
title('Before and After adding Gaussian noise');
xlabel('t');
legend('Xt','Yt=X(t)+W(t)');

%% TASK 8
%Create ÕI(t) and ÕQ(t) signals
Yi_1 = Yt.*(cos(2*pi*Fo*transpose(t_Xtq)));
Yq_1 = Yt.*(-sin(2*pi*Fo*transpose(t_Xtq)));

Tyi1_total = length(t_Xtq)*Ts;
Tyq1_total = length(t_Xtq)*Ts;
%Fourier Transform and compute Periodogram for both signals
Fyi1_x = fftshift(fft(Yi_1,Nf)*Ts);
Pyi1_x = (abs(Fyi1_x).^2)/Tyi1_total;

Fyq1_x = fftshift(fft(Yq_1,Nf)*Ts);
Pyq1_x = (abs(Fyq1_x).^2)/Tyq1_total;
%Show results for each signal
figure();
subplot(2,1,1);
plot(t_Xtq,Yi_1);
title('YI(t)');
xlabel('t');
subplot(2,1,2);
plot(F,Pyi1_x);
title('Periodogram for YI(t)'); 
ylabel('Py(F)');
xlabel('Frequency (Hz)');

figure();
subplot(2,1,1);
plot(t_Xtq,Yq_1);
title('YQ(t)')
xlabel('t');
subplot(2,1,2);
plot(F,Pyq1_x);
title('Periodogram for YQ(t)'); 
ylabel('Py(F)');
xlabel('Frequency (Hz)');

%% TASK 9
%Using SRRC pulse (ph) from Task 3
%Compute Convolution and create time axis.
Yi_2 = conv(ph,Yi_1)*Ts;
t_Yi_2 = (t(1)+t_Xtq(1):Ts:t(end)+t_Xtq(end));
Yq_2 = conv(ph,Yq_1)*Ts;
t_Yq_2 = (t(1)+t_Xtq(1):Ts:t(end)+t_Xtq(end));

%Tail cutting
counter = 0;
for n = 1:length(t_Yi_2)
    if t_Yi_2(n)<0
       counter = counter+1;
    end
end
Yi_cut = Yi_2(counter+1:(length(t_Yi_2)-(counter)));
t_Yi_cut = t_Yi_2(counter+1:(length(t_Yi_2)-(counter)));

Yq_cut = Yq_2(counter+1:(length(t_Yi_2)-(counter)));
t_Yq_cut = t_Yq_2(counter+1:(length(t_Yi_2)-(counter)));
%Fourier Transform and compute Periodogram for both cutted signals
Tyi_2_total = length(t_Yi_cut)*Ts;
Fyi_2 = fftshift(fft(Yi_cut,Nf)*Ts);
Pyi_2 = (abs(Fyi_2).^2)/Tyi_2_total;

Tyq_2_total = length(t_Yq_cut)*Ts;
Fyq_2 = fftshift(fft(Yq_cut,Nf)*Ts);
Pyq_2 = (abs(Fyq_2).^2)/Tyq_2_total;
%Show results for each signal
figure();
subplot(2,1,1);
plot(t_Yi_cut,Yi_cut);
title('Y_i(t) filtered');
subplot(2,1,2);
plot(F,Pyi_2);
title('Periodogram of YI(t) after filtering with SRRC pulse'); 
ylabel('Py(F)');
xlabel('Frequency (Hz)');


figure();
subplot(2,1,1);
plot(t_Yq_cut,Yq_cut);
title('YQ(t) after filtering with SRRC pulse')
xlabel('t');
subplot(2,1,2);
plot(F,Pyq_2);
title('Periodogram of YQ(t) after filtering with SRRC pulse'); 
ylabel('Py(F)');
xlabel('Frequency (Hz)');
%% TASK 10
%Downsampling and creating Yi,n and Yq,n
YI = downsample(Yi_cut,over);
YQ = downsample(Yq_cut,over);
%Creating Yn sequence
Yn=[YI YQ];
%Use Scatterplot for showing results
scatterplot(Yn);
%% TASK 11
[est_X, est_bit_seq] = detect_PSK_8(Yn);
%% TASK 12
Num_Of_Symbol_errors = num_of_symbol_errors(est_X,Xn)
%% TASK 13
Num_Of_Bit_errors = num_of_bit_errors(est_bit_seq,b)
