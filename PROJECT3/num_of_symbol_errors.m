function [ num_of_symbol_errors] = num_of_symbol_errors(est_X,X )
%Compares  the symbols at the imput(X) with the estimation of them at the
%output(est_X)
%       est_X:        Output Symbol Sequence Estimation 
%       X:                Input Symbol Sequence
%Returns total number of mismatches found
num_of_symbol_errors=0;
N=length(X);
%Posible values [-1, -sqrt(2)/2 , 0 , sqrt(2)/2 , 1]
%Multiply by 2 and make values integer so we can compare
est_X_toInt=round(2.*est_X);
X_toInt=round(2.*X);
%Now possible values are [-2 ,-1, 0, 1, 2]
for i=1:N
    %There is a mismatch.Increase total errors
    if est_X_toInt(i,1)~=X_toInt(i,1) || est_X_toInt(i,2)~=X_toInt(i,2)
        num_of_symbol_errors=num_of_symbol_errors+1;
    end
end
end


