function [out01,outpm] = getG1()

% Function getG1 will return the G1 stream
% for generating GPS C/A code.
% Guttorm R. Opshaug 04/04/01
% Penny Axelrad added +1/-1 output 8/29/20


lcode = 1023;
register = ones(1,10);
code = zeros(10,1);
code([3 10]) = ones(2,1);

out01 = zeros(1,lcode);
out01(1) = register(10);

for ij = 2: lcode
   temp = register(1:9);
   register(1) = mod(register*code,2);
   register(2:10) = temp;
   out01(ij) = register(10);
end

outpm = 2*out01 - 1;