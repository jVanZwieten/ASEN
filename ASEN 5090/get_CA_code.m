function [out,outpm] = get_CA_code(PRN)

% Function get_CA_code(PRN) will generate one full epoch
% of C/A code data for the given PRN number
% Guttorm R. Opshaug 04/04/01
% Penny Axelrad added +1/-1 output 8/29/20

S = getS(PRN);
G1 = getG1;
G2 = getG2(S(1),S(2));
out = mod(G1+G2,2);
outpm = 2*out - 1;