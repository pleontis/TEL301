%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all;
close all;
%Data
T=0.1;
over = 10;
a=0.5;
A=5;
N=100;
Ts=T/over;
%C.1
figure(1);
b = (sign (randn(N,1)) + 1)/2;
%C.2.a
x = bits_to_2PAM(b);
%C.2.b
Xd = 1/Ts*upsample(x,over);
t_d = 0:Ts:N/over -Ts;
plot(t_d, Xd);
title('Xd(t)');
xlabel('Xd(t)');
ylabel('t_d');

%C.2.c
figure(2);
[ph,t]=srrc_pulse(T,over,A,a);

tX= (t_d(1)+t(1):Ts:t_d(end)+t(end));
X = conv(Xd, ph).*Ts;
plot(tX,X);

title('Convolution of Xd(t) and PH(t)');
xlabel('tX');
ylabel('X');

%C.2.d
tZ = ( tX(1) + t(1):Ts:tX(end)+t(end));
Z = conv(X,ph).*Ts;
figure(3);
plot(tZ,Z);
title('Convolution of X(t) and pH(t)');
xlabel('tZ');
ylabel('Z(kt)');

figure(4);
plot(tZ,Z);
hold on;
stem([0:N-1]*T,x);
title('Compare Z(kT) with Xd');
xlabel('tZ');
ylabel('Z(kt)');
