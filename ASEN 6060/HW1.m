clc; clear; close all;
format long e
addpath('..')

Gtil = CelestialParameters.universalGravitationalConstant;

%% 1a
GMtil_earth = CelestialParameters.gravityParameter_earth;
Mtil_earth = GMtil_earth/Gtil

GMtil_luna = CelestialParameters.gravityParameter_luna;
Mtil_luna = GMtil_luna/Gtil

mStar_earthLuna = CR3BPUtilities.characteristicMass(Mtil_earth, Mtil_luna)
mu_earthLuna = CR3BPUtilities.massRatio(Mtil_earth, Mtil_luna)

lStar_earthLuna = CelestialParameters.semiMajorAxis_luna;

tStar_earthLuna = CR3BPUtilities.characteristicTime(lStar_earthLuna, mStar_earthLuna)

GMtil_sol = CelestialParameters.gravityParameter_sol;
Mtil_sol = GMtil_sol/Gtil

mStar_solEarth = CR3BPUtilities.characteristicMass(Mtil_sol, Mtil_earth)
mu_solEarth = CR3BPUtilities.massRatio(Mtil_sol, Mtil_earth)

lStar_solEarth = CelestialParameters.semiMajorAxis_earth;

tStar_solEarth = CR3BPUtilities.characteristicTime(lStar_solEarth, mStar_solEarth)

%% 2bi
xVec_0i = [0.98; 0; 0; 0; 1.2; 0];
tSpan_i = [0 2];
tolerancePower = -12;
[t_i, X_i] = CR3BPUtilities.integrateCr3bp(xVec_0i, mu_earthLuna, tSpan_i, tolerancePower);

CR3BPUtilities.plotCr3bpTrajectory(X_i, mu_earthLuna, '', 'Luna');
title('trajectory i')

%% 2bii
xVec_0ii = [0.98; 0; 0; 0; 1.7; 0];
tSpan_ii = [0 8];
[t_ii, X_ii] = CR3BPUtilities.integrateCr3bp(xVec_0ii, mu_earthLuna, tSpan_ii, tolerancePower);

CR3BPUtilities.plotCr3bpTrajectory(X_ii, mu_earthLuna, '', 'Luna');
title('trajectory ii')

%% 2biii
xVec_0iii = [0.12; 0; 0; 0; 3.45; 0];
tSpan_iii = [0 25];
[t_iii, X_iii] = CR3BPUtilities.integrateCr3bp(xVec_0iii, mu_earthLuna, tSpan_iii, tolerancePower);

CR3BPUtilities.plotCr3bpTrajectory(X_iii, mu_earthLuna, 'Earth', 'Luna');
title('trajectory iii')

%% 2biv
xVec_0iv = [0.12; 0; 0; 0; 3.48; 0];
tSpan_iv = [0 25];
[t_iv, X_iv] = CR3BPUtilities.integrateCr3bp(xVec_0iv, mu_earthLuna, tSpan_iv, tolerancePower);

CR3BPUtilities.plotCr3bpTrajectory(X_iv, mu_earthLuna, 'Earth', 'Luna');
title('trajectory iv')

%% 2c
C_i = jacobiConstantThroughoutTrajectory(X_i, mu_earthLuna);
C_ii = jacobiConstantThroughoutTrajectory(X_ii, mu_earthLuna);
C_iii = jacobiConstantThroughoutTrajectory(X_iii, mu_earthLuna);
C_iv = jacobiConstantThroughoutTrajectory(X_iv, mu_earthLuna);

dC_i = Utilities.differenceAtEachStep(C_i);
dC_ii = Utilities.differenceAtEachStep(C_ii);
dC_iii = Utilities.differenceAtEachStep(C_iii);
dC_iv = Utilities.differenceAtEachStep(C_iv);

dC_iMax = max(dC_i)
dC_iiMax = max(dC_ii)
dC_iiiMax = max(dC_iii)
dC_ivMax = max(dC_iv)

figure;
tiledlayout(4, 1);

nexttile;
plot(t_i(2:end), dC_i);
title('dC_i vs t_i');
xlabel('Time');
ylabel('dC_i');

nexttile;
plot(t_ii(2:end), dC_ii);
title('dC_{ii} vs t_{ii}');
xlabel('Time');
ylabel('dC_{ii}');

nexttile;
plot(t_iii(2:end), dC_iii);
title('dC_{iii} vs t_{iii}');
xlabel('Time');
ylabel('dC_{iii}');

nexttile;
plot(t_iv(2:end), dC_iv);
title('dC_{iv} vs t_{iv}');
xlabel('Time');
ylabel('dC_{iv}');


%% 2d
vStar_earthLuna = lStar_earthLuna/tStar_earthLuna;
xVec_0iii
xVecTilde_0iii = [lStar_earthLuna*xVec_0iii(1:3); vStar_earthLuna*xVec_0iii(4:6)]
tSpanTilde_iii = tSpan_iii*tStar_earthLuna
tSpanTilde_iiiDays = Utilities.s2days(tSpanTilde_iii)

r_1 = lStar_earthLuna*norm(xVec_0iii(1:3) - [-mu_earthLuna; 0; 0])
tSpan_iiiPeriods = tSpan_iii/2/pi

%% 3b
[t_3b, X_3b, t_f3b, xVec_f3b, index_f3b] = CR3BPUtilities.integrateCr3bpUntil(xVec_0iii, mu_earthLuna, tSpan_iii, tolerancePower, @crossXAxis);
xVec_f3b
t_f3b

%% 4a
CR3BPUtilities.plotZeroVelocityCurve(3.189, mu_earthLuna, 'Earth', 'Luna'); title("ZVC C = 3.189");
CR3BPUtilities.plotZeroVelocityCurve(3.173, mu_earthLuna, 'Earth', 'Luna'); title("ZVC C = 3.173");
CR3BPUtilities.plotZeroVelocityCurve(3.013, mu_earthLuna, 'Earth', 'Luna'); title("ZVC C = 3.013");
CR3BPUtilities.plotZeroVelocityCurve(2.995, mu_earthLuna, 'Earth', 'Luna'); title("ZVC C = 2.995");

%% 4d
x = .84;
C = CR3BPUtilities.jacobiConstant([x; 0; 0], mu_earthLuna);
CR3BPUtilities.plotZeroVelocityCurve(C, mu_earthLuna, 'Earth', 'Luna');
title(['ZVC x = ' num2str(x) ' C = ' num2str(C, 3)]);
xlim([0.8 0.88]); ylim([-0.04 0.04]);

x = 1.16;
C = CR3BPUtilities.jacobiConstant([x; 0; 0], mu_earthLuna);
CR3BPUtilities.plotZeroVelocityCurve(C, mu_earthLuna, 'Earth', 'Luna');
title(['ZVC x = ' num2str(x) ' C = ' num2str(C, 3)]);
xlim([1.1 1.2]); ylim([-0.05 0.05]);

x = .5;
y = .83:.01:.88;
for y_i = y
    C = CR3BPUtilities.jacobiConstant([x; y_i; 0], mu_earthLuna);
    CR3BPUtilities.plotZeroVelocityCurve(C, mu_earthLuna, 'Earth', 'Luna');
    title(['ZVC r =[' num2str(x) '; ' num2str(y_i) '] C = ' num2str(C, 3)]);
end

function [value,isterminal,direction] = crossXAxis(t, XVec)
    value = XVec(2);
    isterminal = 1;
    direction = 1;
end

function CVec = jacobiConstantThroughoutTrajectory(XMat, mu)
    CVec = zeros(size(XMat, 1), 1);
    for i = 1:size(XMat, 1)
        CVec(i) = CR3BPUtilities.jacobiConstant(XMat(i, :)', mu);
    end
end