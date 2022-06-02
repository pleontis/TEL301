function [ x ] = bits_to_PSK_8( bit_seq )
%bits_to_PSK_8 converts a bit sequence into 8-PSK coding
%Picked Gray Code:
%               000->001->011->010->110->100->101->111
    for i=1:length(bit_seq)
        if bit_seq(i,1)==0 && bit_seq(i,2)==0 && bit_seq(i,3)==0
            x(i,1)=cos(2*pi*0/8);
            x(i,2)=sin(2*pi*0/8);
        elseif bit_seq(i,1)==0 && bit_seq(i,2)==0 && bit_seq(i,3)==1
            x(i,1)=cos(2*pi*1/8);
            x(i,2)=sin(2*pi*1/8);
         elseif bit_seq(i,1)==0 && bit_seq(i,2)==1 && bit_seq(i,3)==1
            x(i,1)=cos(2*pi*2/8);
            x(i,2)=sin(2*pi*2/8);
         elseif bit_seq(i,1)==0 && bit_seq(i,2)==1 && bit_seq(i,3)==0
            x(i,1)=cos(2*pi*3/8);
            x(i,2)=sin(2*pi*3/8);
         elseif bit_seq(i,1)==1 && bit_seq(i,2)==1 && bit_seq(i,3)==0
            x(i,1)=cos(2*pi*4/8);
            x(i,2)=sin(2*pi*4/8);
         elseif bit_seq(i,1)==1 && bit_seq(i,2)==0 && bit_seq(i,3)==0
            x(i,1)=cos(2*pi*5/8);
            x(i,2)=sin(2*pi*5/8);
         elseif bit_seq(i,1)==1 && bit_seq(i,2)==0 && bit_seq(i,3)==1
            x(i,1)=cos(2*pi*6/8);
            x(i,2)=sin(2*pi*6/8);
        elseif bit_seq(i,1)==1 && bit_seq(i,2)==1 && bit_seq(i,3)==1
            x(i,1)=cos(2*pi*7/8);
            x(i,2)=sin(2*pi*7/8);
        end
    end
end

