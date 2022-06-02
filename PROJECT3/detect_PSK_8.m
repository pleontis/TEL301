function [est_X,exit_bit_seq ] = detect_PSK_8( Y )
%Compare distances and decide according to nearest neighbor rule
J=sqrt(2)/2;
for i=1:length(Y)
        dx0 = norm([1 0]-Y(i,:));
        dx1 = norm([J J]-Y(i,:));
        dx2 = norm([0 1]-Y(i,:));
        dx3 = norm([-J J]-Y(i,:));
        dx4 = norm([-1 0]-Y(i,:));
        dx5 = norm([-J -J]-Y(i,:));
        dx6 = norm([0 -1]-Y(i,:));
        dx7 = norm([J -J]-Y(i,:));

        dmin = min([dx0,dx1,dx2,dx3,dx4,dx5,dx6,dx7]);
% Decision of symbol m=0,1,..7 
        if (dmin == dx0)
            est_X(i,:) = [1 0];                
        elseif (dmin == dx1)              
            est_X(i,:) = [J J];                 
        elseif (dmin == dx2)
            est_X(i,:) = [0 1];                 
        elseif (dmin == dx3)
            est_X(i,:) = [-J J];
        elseif (dmin == dx4)
            est_X(i,:) = [-1 0];
        elseif (dmin == dx5)
            est_X(i,:) = [-J -J];
        elseif (dmin == dx6)
            est_X(i,:) = [0 -1];
        elseif (dmin == dx7)
            est_X(i,:) = [J -J];
        end
end

%Decoding
%Gray code picked for bits_to_PSK_8 :
%000->001->011->010->110->100->101->111
for i=1:length(Y)
    if est_X(i,1)==1  && est_X(i,2)==0                                        %m=0
        exit_bit_seq(i,:)=[0 0 0];
    elseif est_X(i,1)==J  && est_X(i,2)==J                                  %m=1
        exit_bit_seq(i,:)=[0 0 1];
    elseif est_X(i,1)==0  && est_X(i,2)==1                                 %m=2
        exit_bit_seq(i,:)=[0 1 1 ];
    elseif est_X(i,1)==-J  && est_X(i,2)==J                                 %m=3
        exit_bit_seq(i,:)= [ 0 1 0];
    elseif est_X(i,1)==-1  && est_X(i,2)==0                                %m=4
        exit_bit_seq(i,:)=[1 1 0];
    elseif est_X(i,1)==-J && est_X(i,2)==-J                                 %m=5
        exit_bit_seq(i,:)=[1 0 0];
    elseif est_X(i,1)==0  && est_X(i,2)==-1                                %m=6
        exit_bit_seq(i,:)=[1 0 1];
    elseif est_X(i,1)==J  && est_X(i,2)==-J                                 %m=7                                                
        exit_bit_seq(i,:)=[1 1 1];
    end
end
end