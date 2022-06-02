function [ x] = bits_to_4PAM( b1, b2 )
    for i=1:1:length(b1)
        if(b1(i)==0)
             if(b2(i)==0)
                x(i)=+3;
             else
                x(i)=+1;
             end
        else %That means that b1==1
            if(b2(i)==0)
                x(i)=-3;
            else
                x(i)=-1;
            end
        end
    end
end

