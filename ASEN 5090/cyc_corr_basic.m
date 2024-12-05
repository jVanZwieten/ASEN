function [Rxy,lag] = cyc_corr_basic(x,y)
% Function to compute the cyclic correlation for ASEN 5090
% First input is the "received" signal
% Second input is the "replica" signal - this is the one that shifts

n = length(x);
x = reshape(x,1,n);         % make sure x is a row vector
yshifted = reshape(y,n,1);  % make the replica of y a column vector

lag = [0:n]; Rxy = NaN(size(lag));

for i = 1:length(lag)  % sweep the lag from 0 to n
    Rxy(i) = x*yshifted;
    yshifted = [yshifted(end); yshifted(1:end-1)];
end
