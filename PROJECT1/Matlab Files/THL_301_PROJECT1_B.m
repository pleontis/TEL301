%Polydoros Prinitis -Panagiotis Leontis
%2018030098 - 2018030099
clear all;
close all;
%Data
T=10^(-3);
over=10;
Ts=T/over;
A=5;

rolloff = [0 0.5 1];
for k=[0 1 2 3] %For each k
    
    for i=[1 2 3]   %For each a.
        
         %b1
         figure;    %New figure for a set of k and a
         subplot(2,1,1);    %sublpot for G(t) & G(t-kT)
         [ph1,t] = srrc_pulse(T,over,A,rolloff(i));     
         %Create the G(t-kT);
         pha= [ph1, zeros(1,10*k)];
         phb=[zeros(1,10*k),ph1];
         %Create the new time after adding the zeros
         ts=[0:Ts:length(t)*Ts+T*k-Ts];
         %Plot G(t)
         plot1 = plot(ts,pha);
         hold on;
         %Plot G(t-kT)
         plot2 = plot(ts,phb);
         
         title(['PH(t) & PH(t-kT) Pulses  a= ' ,num2str(rolloff(i)), '  k= ',num2str(k)]);
         ylabel('PH(t) & PH(t-kT) ');
         xlabel('ts');
         legend([plot1, plot2],'PH(t)' , 'PH(t-kT)');
         
         %b2
         subplot(2,1,2);    %subplot for G(t) * G(t-kT)
         phR = pha.*phb;
         plot(ts, phR);
         
         title(['PH(t) * PH(t-kT)   a= ' ,num2str(rolloff(i)), '  k= ',num2str(k)]);
         ylabel('PH(t) * PH(t-kT) ');
         xlabel('ts');
         
         %b3
         integral = sum(phR).*Ts;   %Calculate integral
         fprintf('\nIntegral of PH(t)*PH(t-kT) for a= %2f and k=%d is : %2f  ',rolloff(i),k,integral);
         
    end
end