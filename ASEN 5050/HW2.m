format shortG
% 1
clear
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

Omega = astroUtilities.RAANFromN(N)
fprintf('%.6g\n', rad2deg(Omega))

omega = -astroUtilities.argumentOfPeriapsisFromNE(N, E)
fprintf('%.6g\n', rad2deg(omega))

nu = -astroUtilities.trueAnomalyFromRE(R, E)
fprintf('%.6g\n', rad2deg(nu))

theta = omega + nu
fprintf('%.6g\n', rad2deg(theta))

% b
C = astroUtilities.DirectionCosineMatrix(Omega, omega, nu, i)
Ct = C'
Ct*C
C*Ct

R_radialFrame = Ct*R
V_radialFrame = Ct*V

phi = astroUtilities.flightPathAngleFromENu(e, nu)
fprintf('%.6g\n', rad2deg(phi))

V_radialFrame_0 = [v*sin(phi); v*cos(phi); 0]

% c
r_p = astroUtilities.periapsisFromAE(a, e)
R_p = astroUtilities.PeriapsisFromRE(r_p, E)

v_p = astroUtilities.velocityAtApsisFromHR(h, r_p)
astroUtilities.specificMechanicalNRGFromrv(r_p, v_p, mu_mars)

C_p = astroUtilities.DirectionCosineMatrix(Omega, omega, 0, i)

V_p = C_p*[0; v_p; 0]
cross(R_p, V_p)
norm(cross(R_p, V_p))

% % 2
% clear
% mu_moon = 4.9902799e3; % km^3/s^2
% r_moonEquator = 1738; % km

% Nhat = [0.6428; -0.766; 0];
% Hhat = [-0.3237; -0.2717; 0.9063];
% E = [0.0475; 0.3755; 0.1295];

% % a
% e = norm(E)
% i = astroUtilities.inclinationFromHhat(Hhat)
% fprintf('%.6g\n', rad2deg(i))
% Omega = -astroUtilities.RAANFromNhat(Nhat)
% fprintf('%.6g\n', rad2deg(Omega))
% omega = astroUtilities.argumentOfPeriapsisFromNE(Nhat, E)
% fprintf('%.6g\n', rad2deg(omega))

% Ehat = E/e;
% omega2 = acos(dot(Nhat, Ehat))

% % b
% r = 4070.6; % km

% a = r*(1 + e*cos(-omega))/(1 - e^2)