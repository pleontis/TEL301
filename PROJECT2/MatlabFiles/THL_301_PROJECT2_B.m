%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all; 
close all;
%% B4
T=10^(-2);
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=1/2;
Ns=2048;

 F=-Fs/2:Fs/Ns:Fs/2-Fs/Ns;
%Choosing a random value for fo in space 50<fo<450
fo=100;
%Create pulse
[ph,t] = srrc_pulse(T,over,A,a);
phF=fftshift(fft(ph,Ns).*Ts);
%Experiment for k repetitions
N=100;
t_xN=(0:Ts:N/(1/T)-Ts);
t_conv=t(1)+t_xN(1):Ts:t(end)+t_xN(end);
len_t=length(t_conv)*Ts; 
%Change direction of array so we can multiply
to=t_conv.';

k=500;
for i=1:k
    %Create 2-PAM sequence
    bs = (sign (randn(N,1)) + 1)/2;
    x = bits_to_2PAM(bs);
    xN=(1/Ts)*upsample(x,over);
    %Compute Convolution and create time axis
    x_conv=conv(ph,xN).*Ts;
    %Set random variable theta with uniform distribution and create Y(t) signal
    theta=unifrnd(0,2*pi);
    yt=x_conv.*cos(2.*pi.*fo.*to+theta);
    %Compute Fourier Transform and Power Spectral Density
    yF=fftshift(fft(yt,Ns)).*Ts;
    Py=(abs(yF).^2)/len_t;
    P(i,:)=Py;
end
Pavg=sum(P)/k;
%Semilogy results
semilogy(F,Pavg);
title(['Experimental Diagram for k= ',num2str(k)]);
xlabel('Frequency');
ylabel('Px(F)')



