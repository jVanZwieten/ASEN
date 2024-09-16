clear

% 1
mu_mars = 4.305e4; % km^3/s^2
r_marsEquator = 3397.2; % km
R = [3.62067e3;
-3.19925e2;
-4.20645e2]; % km
V = [-4.28843e-1;
-3.00176e-2;
-3.39801]; % km/s

% a
r = norm(R)
v = norm(V)

H = astroUtilities.angularMomentumFromRV(R, V);
h = norm(H)

specMechNRG = astroUtilities.specificMechanicalNRGFromRV(R, V, mu_mars)
a = astroUtilities.semiMajorAxisFromSpecificMechanicalNRG(specMechNRG, mu_mars)

E = astroUtilities.EccentricityFromRV(R, V, mu_mars)
e = norm(E)

i = astroUtilities.inclinationFromH(H)
fprintf('%.6g', rad2deg(i))

N = astroUtilities.LineOfNodesFromH(H)
n = norm(N)

Omega = astroUtilities.longitudeOfAscendingNodeFromN(N)
fprintf('%.6g\n', rad2deg(Omega))

omega = astroUtilities.argumentOfPeriapsisFromNE(N, E);
omega = 2*pi - omega
fprintf('%.6g\n', rad2deg(omega))

nu = astroUtilities.trueAnomalyFromRE(R, E);
nu = 2*pi - nu
fprintf('%.6g\n', rad2deg(nu))
