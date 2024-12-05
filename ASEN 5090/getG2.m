function out21 = getG2(S1,S2)

% Function getG2(S1,S2) will return the G21 stream
% for generating the GPS C/A code.
% S1 and S2 are tap positions for the G2 shift register
% Guttorm R. Opshaug 04/04/01

lcode = 1023;
register = ones(1,10);
temp = zeros(1,9);
code = zeros(10,1);
code([2 3 6 8 9 10]) = ones(6,1);

out21 = zeros(1,lcode);
out21(1) = mod(register(S1) + register(S2),2);

for ij = 2: lcode
   temp = register(1:9);
   register(1) = mod(register*code,2);
   register(2:10) = temp;
   out21(ij) = mod(register(S1) + register(S2),2);
end