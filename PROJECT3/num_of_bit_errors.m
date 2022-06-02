function [ num_of_bit_errors ] = num_of_bit_errors(est_bit_seq,b )
%Compares the bit sequence of input signal and the estimated bit sequence of output
%       est_bit_seq:        Output Binary 3-bit Sequence 
%       b:                           Input Binary 3-bit Sequence
%Returns total number of mismatches found

%Number or triads
N=length(b);
num_of_bit_errors=0;
%Repeat for all triads
for i=1:N
    %There is a mismatch. Increase total errors
    if b(i,1)~=est_bit_seq(i,1) || b(i,2)~=est_bit_seq(i,2) || b(i,3)~=est_bit_seq(i,3)
        num_of_bit_errors=num_of_bit_errors+1;
    end
end

