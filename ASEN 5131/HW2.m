clear
close all

%% 1b
gamma = 1.4;
M = 0:10;
KEh0 = HypersonicsUtilities.kineticEnergyRatio(M, gamma);

function KEh0 = kineticEnergyRatio(M, gamma)
    KEh0 = M.^2./(2/(gamma - 1) + M.^2);
end

figure
plot(M, KEh0)
xlabel("M")
ylabel("Kinetic to Total Energy Ratio")
title("Kinetic to Total Energy Ratio Profile for \gamma=1.4")

%% 4a
h_5000 = 9866.70; % kJ/kg
h_8000 = 38039.6; % kJ/kg
h_ratio = h_8000/h_5000
h_ratioWOGasEffects = 8000/5000

%% 4d
h_5000 = 7867.99;
h_8000 = 18171.3;
h_ratio = h_8000/h_5000

N2 = [7.5203e-1 1.1340e-1 7.2562e-1 6.2901e-1];
O2 = [2.6686e-3 1.8246e-5 6.5880e-2 2.1242e-3];
N2Ratio_1atm = N2(2)/N2(1)
N2Ratio_100atm = N2(4)/N2(3)
O2Ratio_1atm = O2(2)/O2(1)
O2Ratio_100atm = O2(4)/O2(3)