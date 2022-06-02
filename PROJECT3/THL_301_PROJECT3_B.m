%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all;
close all;

N=100;
SNR=[-2:2:16];
T=10^(-3);
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=1/2;
Fo=2000;
K=200;
%Create SRRC pulse 
[ph,t] = srrc_pulse(T,over,A,a);

%For different values of SNR repeat experiment and calculate
%Symbol Error Probability and Bit Error Probability
for test=1:length(SNR)
    %Initalize for each repetition
    Num_Of_Symbol_Errors=0;
    Num_Of_Bits_Errors=0;
    
    for samples=1:K
        b = (sign (randn(N,3)) + 1)/2;
        
        Xn=bits_to_PSK_8(b);
        XI=Xn(:,1);
        XQ=Xn(:,2);
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
 
        %Create XI(t) and XQ(t) signals (modulation)
        XTi = 2*(Xti).*(cos(2*pi*Fo*transpose(t_Xti)));

        XTq = -2*(Xtq).*(sin(2*pi*Fo*transpose(t_Xtq)));
        
       %Create channel input signal
        Xt = XTi+XTq;
        %Compute variance and create Gaussian Noise
        s2w = 1/(Ts*(10^(SNR(test)/10)));
        s2n=Ts*s2w/2;
        Wt = sqrt(s2w).*randn(length(Xt),1);
        %Create Y(t) signal
        Yt = Xt+Wt;
        
        %Create ÕI(t) and ÕQ(t) signals( de-modulation)
        Yi_1 = Yt.*(cos(2*pi*Fo*transpose(t_Xtq)));
        Yq_1 = Yt.*(-sin(2*pi*Fo*transpose(t_Xtq)));
        
        %Compute Convolution and create time axis.
        Yi_2 = conv(ph,Yi_1)*Ts;
        t_Yi_2 = (t(1)+t_Xtq(1):Ts:t(end)+t_Xtq(end));
        Yq_2 = conv(ph,Yq_1)*Ts;

        %Tail cutting
        counter = 0;
        for n = 1:length(t_Yi_2)
          if t_Yi_2(n)<0
            counter = counter+1;
          end
        end
        Yi_cut = Yi_2(counter+1:(length(t_Yi_2)-(counter)));
        Yq_cut = Yq_2(counter+1:(length(t_Yi_2)-(counter)));
        
       %Downsampling and creating Yi,n and Yq,n
        YI = downsample(Yi_cut,over);
        YQ = downsample(Yq_cut,over);
        %Creating Yn sequence
        Yn=[YI YQ];

        %Detect 8-PSK Sequence and calculate symbol errors and bit errors
        [est_X,est_bit_seq]=detect_PSK_8(Yn);
        
        Num_Of_Symbol_Errors=Num_Of_Symbol_Errors+num_of_symbol_errors(est_X,Xn);
        Num_Of_Bits_Errors= Num_Of_Bits_Errors+num_of_bit_errors(est_bit_seq,b);
    end
    %Store probability of error on each sequence into matrix for every repetition
    Symbol_Errors(1,test)=Num_Of_Symbol_Errors/(N*K);
    Bits_Errors(1,test)=Num_Of_Bits_Errors/(N*K*3);
    
    symbol_bound(test)=2*Q(1./(sqrt(s2n))*sin(pi/8));
    %Log2(8)=3
    bit_bound(test)=symbol_bound(test)/3;
end

%Show results with semilogy
figure()
semilogy(SNR,symbol_bound,'red');
hold on;
semilogy(SNR,Symbol_Errors,'blue');
hold off;
xlabel('SNR (db)');
ylabel('Symbol Error Probability');
legend('Smart Upper Bound','Estimated Symbol Error Probability'); 
title('Monte Carlo Method for SER approximation');

figure()
semilogy(SNR,bit_bound,'red');
hold on;
semilogy(SNR,Bits_Errors,'blue');
hold off;
xlabel('SNR (db)');
ylabel('Bit Error Probability');
legend(' Bound','Estimated Bit Error Probability');
title('Monte Carlo Method for BER approximation');