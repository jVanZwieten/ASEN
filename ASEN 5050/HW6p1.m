addpath("..")
clear
close all

AU = astroUtilities;
LU = lambertUtilities;

mu_sun = CelestialParameters.gravityParameter_sun;

t_1 = 2457754.5; % Julian Date (days)
t_2 = 2457871.5; % Julian Date (days)

rVec_1 = [-2.686982e7; 1.326980e8; 5.752566e7]; % km
vVec_1 = [-29.781722; -5.101774; -2.210394]; % km/s
rVec_2 = [-5.642901e7; -8.571048e7; -3.499466e7]; % km
vVec_2 = [29.658341; -16.0911; -9.116674]; % km/s

longWay = false;


TOF = t_2 - t_1
TOF = Utilities.days2s(TOF)

[r_n, rHat_n] = Utilities.magnitudeDirection([rVec_1 rVec_2]);
r_1 = r_n(1)
r_2 = r_n(2)
rHat_1 = rHat_n(:, 1)
rHat_2 = rHat_n(:, 2)

deltaNu = Utilities.normalizeAngle(Utilities.angleBetweenUnitVectors(rHat_1, rHat_2))

[c, s] = LU.lambertGeometricQuantities(r_1, r_2, deltaNu)

TOF_parabolic = LU.lambertsEquationParabolic(c, s, mu_sun, longWay)
hyperbolicTransfer = TOF < TOF_parabolic

TOF_min = LU.lambertsEquationMin(c, s, longWay, mu_sun)
longTOF = TOF > TOF_min

toleranceFactor = 10^-10;
a_guess = s/2 + 10^6
a = lambertIterator(c, s, mu_sun, longTOF, longWay, hyperbolicTransfer, TOF, a_guess, toleranceFactor)

[~, alpha, beta] = LU.lambertsEquation(c, s, a, mu_sun, longTOF, longWay);
p = LU.semiLatusRectum(a, c, s, r_1, r_2, alpha, beta)
e = AU.eccentricityFromSemiLatusRectum(p, a)

nu_1 = AU.trueAnomalyFrompre(p, r_1, e)
nu_2 = AU.trueAnomalyFrompre(p, r_2, e)

nu_1 = -nu_1;
nu_2 = -nu_2;

f = AU.fFunction(r_2, p, deltaNu)
g = AU.gFunction(r_2, r_1, mu_sun, p, deltaNu)
vVec_x1 = (rVec_2 - f*rVec_1)/g

fDot = AU.fDotFunction(mu_sun, p, deltaNu, r_2, r_1)
gDot = AU.gDotFunction(r_1, p, deltaNu)
vVec_x2 = fDot*rVec_1 + gDot*vVec_x1

deltaVVec_1 = vVec_x1 - vVec_1
deltaV_1 = norm(deltaVVec_1)
deltaVVec_2 = vVec_2 - vVec_x2
deltaV_2 = norm(deltaVVec_2)

% [a, e, nu_1, nu_2, vVec_x1, vVec_x2] = lambertSolver(rVec_1, rVec_2, TOF, longWay, mu_sun)