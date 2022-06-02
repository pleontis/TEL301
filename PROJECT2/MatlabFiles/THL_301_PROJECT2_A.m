%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
 clear all;
 close all;
%% A1
T=10^(-2);
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=1/2;
Ns=4096;
F=-Fs/2:Fs/Ns:Fs/2-Fs/Ns;
[ph,t] = srrc_pulse(T,over,A,a);

phF=fftshift(fft(ph,Ns).*Ts);

figure()
 p= semilogy(F,abs(phF).^2);
 title('Semilogy for Power Spectral Density');
xlabel('F');
ylabel('|PH(F)|^2');

%% A2
N=100;
bs = (sign (randn(N,1)) + 1)/2;
x = bits_to_2PAM(bs);
xN=(1/Ts)*upsample(x,over);
t_xN=(0:Ts:N/(1/T)-Ts);
%Compute Convolution and create time axis
x_conv=conv(ph,xN).*Ts;
t_conv=t(1)+t_xN(1):Ts:t(end)+t_xN(end);

figure();
plot(t_conv,x_conv);
xlabel('t');
ylabel('X(t)');
title('Convolution of Xn and ph(t) 2-PAM');

%% A3 Part 1
len_total=length(t_conv)*Ts;
%Fourier Transform and compute Periodogram
Fx=fftshift(fft(x_conv,Ns)*Ts);
Px=(abs(Fx).^2)/len_total;
%Show results with semilogy
figure()
semilogy(F,Px);
xlabel('Frequency (Hz)');
ylabel('Px(F)');
title('Px(F) with semilogy');
%Show results with plot
figure()
plot(F,Px);
xlabel('Frequency (Hz)');
ylabel('Px(F)');
title('Px(F) with plot');
%% A3 Part 2
k=500;
%Experiment for k repetitions
for i=1:k
    b = (sign(randn(N,1))+1)/2;
    x_test = bits_to_2PAM(b);
    xn=(1/Ts)*upsample(x_test,over);
    x_conv_test=conv(ph,xn)*Ts;
    t_conv_test=t(1)+t_xN(1):Ts:t(end)+t_xN(end);
    len_t=length(t_conv_test)*Ts;
    Fx_test=fftshift(fft(x_conv_test,Ns)*Ts);
    Px=(abs(Fx_test).^2)/len_t;
    P(i,:)=Px;
end
%Comute experimet's average result
Pavg=sum(P)/k;
%Compute theoretical Power Spectral Density
Sx=(var(x)).*abs((phF).^2)./T;
%Show theoretical and experimental results on same plot
figure()
semilogy(F,Sx,'red');
xlabel('Frequency (Hz)');
ylabel('Px(F) and Sx(f)');
title(['Experimental and Theoretical Diagram for k= ',num2str(k)]);
grid on;
hold on;
semilogy(F,Pavg,'blue');
hold off;

%% A4 Part 1
%Create sequence 4-PAM
 N4=N/2;
 b1s=(sign(randn(N4,1))+1)/2;
 b2s=(sign(randn(N4,1))+1)/2;
 x4=bits_to_4PAM(b1s,b2s);
 xN4=(1/Ts)*upsample(x4,over);
 t_xN4=(0:Ts:N4/(1/T)-Ts);
 %Compute Convolution and create time axis
 x_conv=conv(ph,xN4).*Ts;
 t_conv=t(1)+t_xN4(1):Ts:t(end)+t_xN4(end);
 %Show results
 figure();
 plot(t_conv,x_conv);
 xlabel('t');
 ylabel('X(t)')
 title('Convolution of Xn and ph(t) 4-PAM');
 %% A4 Part 2
 k=500;
 %Experiment for k repetitions
 for i=1:k
     b1=(sign(randn(N4,1))+1)/2;
     b2=(sign(randn(N4,1))+1)/2;
     x4_test=bits_to_4PAM(b1,b2);
     xn4=(1/Ts)*upsample(x4_test,over);
     x_conv_test4=conv(ph,xn4)*Ts;
     t_conv_test4=t(1)+t_xN4(1):Ts:t(end)+t_xN4(end);
     len_t4=length(t_conv_test4)*Ts;
     Fx4_test=fftshift(fft(x_conv_test4,Ns)*Ts);
     Px=(abs(Fx4_test).^2)/len_t4;
     P(i,:)=Px;
 end
 %Comute experimet's average result
 Pavg=sum(P)/k;
 %Compute theoretical Power Spectral Density
 Sx=(var(x4)).*abs((phF).^2)./T;
 %Show theoretical and experimental results on same plot
 figure()
semilogy(F,Sx,'red');
xlabel('Frequency (Hz)');
ylabel('Px(F) and Sx(f)');
title(['Experimental and Theoretical Diagram for k= ',num2str(k)]);
grid on;
hold on;
semilogy(F,Pavg,'blue');
hold off;
%% A5
N=100;
over_5=2*over;
T_5=2*T;
Ns=4096;
F=-Fs/2:Fs/Ns:Fs/2-Fs/Ns;
[ph,t] = srrc_pulse(T_5,over_5,A,a);
phF=fftshift(fft(ph,Ns).*Ts);
xN=(1/Ts)*upsample(x,over_5);
t_xN=(0:Ts:N/(1/T_5)-Ts);
%Compute Convolution and create time axis
x_conv=conv(ph,xN).*Ts;
t_conv=t(1)+t_xN(1):Ts:t(end)+t_xN(end);

len_total=length(t_conv)*Ts;
%Fourier Transform and compute Periodogram
Fx=fftshift(fft(x_conv,Ns)*Ts);
Px=(abs(Fx).^2)/len_total;
%Show results with semilogy
figure()
semilogy(F,Px);
xlabel('Frequency (Hz)');
ylabel('Px(F)');
title('Px(F) with semilogy with T"=2T ');
%Show results with plot
figure()
plot(F,Px);
xlabel('Frequency (Hz)');
ylabel('Px(F)');
title('Px(F) with plot T"=2T');
k=500;
%Experiment for k repetitions
for i=1:k
    b = (sign(randn(N,1))+1)/2;
    x_test = bits_to_2PAM(b);
    xn=(1/Ts)*upsample(x_test,over_5);
    x_conv_test=conv(ph,xn)*Ts;
    Fx_test=fftshift(fft(x_conv_test,Ns)*Ts);
    Px=(abs(Fx_test).^2)/len_total;
    P(i,:)=Px;
end
%Comute experimet's average result
Pavg=sum(P)/k;
%Compute theoretical Power Spectral Density
Sx=(var(x)).*abs((phF).^2)./T_5;
%Show theoretical and experimental results on same plot
figure()
semilogy(F,Sx,'red');
xlabel('Frequency (Hz)');
ylabel('Px(F) and Sx(f)');
title(['Experimental and Theoretical Diagram for k= ',num2str(k)]);
grid on;
hold on;
semilogy(F,Pavg,'blue');
hold off;

