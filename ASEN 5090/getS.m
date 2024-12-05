function S = getS(PRN);

% Function getS will return a vector, S, containing
% the 2 taps for generating the G21 stream.
% Returns vector containing S1, S2 for a given PRN number.
% Guttorm R. Opshaug 04/04/01

if((PRN < 1)|(PRN > 37))
   error('PRN numbers between 1 and 37 only!');
end
S = zeros(1,2);

if (PRN == 1)
   S(1) = 2; S(2) = 6;
elseif (PRN == 2)
   S(1) = 3; S(2) = 7;
elseif (PRN == 3)
   S(1) = 4; S(2) = 8;
elseif (PRN == 4)
   S(1) = 5; S(2) = 9;
elseif (PRN == 5)
   S(1) = 1; S(2) = 9;
elseif (PRN == 6)
   S(1) = 2; S(2) = 10;
elseif (PRN == 7)
   S(1) = 1; S(2) = 8;
elseif (PRN == 8)
   S(1) = 2; S(2) = 9;
elseif (PRN == 9)
   S(1) = 3; S(2) = 10;
elseif (PRN == 10)
   S(1) = 2; S(2) = 3;
elseif (PRN == 11)
   S(1) = 3; S(2) = 4;
elseif (PRN == 12)
   S(1) = 5; S(2) = 6;
elseif (PRN == 13)
   S(1) = 6; S(2) = 7;
elseif (PRN == 14)
   S(1) = 7; S(2) = 8;
elseif (PRN == 15)
   S(1) = 8; S(2) = 9;
elseif (PRN == 16)
   S(1) = 9; S(2) = 10;
elseif (PRN == 17)
   S(1) = 1; S(2) = 4;
elseif (PRN == 18)
   S(1) = 2; S(2) = 5;
elseif (PRN == 19)
   S(1) = 3; S(2) = 6;
elseif (PRN == 20)
   S(1) = 4; S(2) = 7;
elseif (PRN == 21)
   S(1) = 5; S(2) = 8;
elseif (PRN == 22)
   S(1) = 6; S(2) = 9;
elseif (PRN == 23)
   S(1) = 1; S(2) = 3;
elseif (PRN == 24)
   S(1) = 4; S(2) = 6;
elseif (PRN == 25)
   S(1) = 5; S(2) = 7;
elseif (PRN == 26)
   S(1) = 6; S(2) = 8;
elseif (PRN == 27)
   S(1) = 7; S(2) = 9;
elseif (PRN == 28)
   S(1) = 8; S(2) = 10;
elseif (PRN == 29)
   S(1) = 1; S(2) = 6;
elseif (PRN == 30)
   S(1) = 2; S(2) = 7;
elseif (PRN == 31)
   S(1) = 3; S(2) = 8;
elseif (PRN == 32)
   S(1) = 4; S(2) = 9;
elseif (PRN == 33)
   S(1) = 5; S(2) = 10;
elseif (PRN == 34)
   S(1) = 4; S(2) = 10;
elseif (PRN == 35)
   S(1) = 1; S(2) = 7;
elseif (PRN == 36)
   S(1) = 2; S(2) = 8;
elseif (PRN == 37)
   S(1) = 4; S(2) = 10;
end