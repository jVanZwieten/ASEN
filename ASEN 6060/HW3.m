clc; clear; close all;
addpath('..')

mu_earthLuna = CR3BPUtilities.mu_earthLuna;
lagrangePoints_earthLuna = CR3BPUtilities.lagrangePoints(mu_earthLuna);
rVec_L1 = lagrangePoints_earthLuna{1};
dRVec_0 = [-.001; 0];
lagrangeStabilities_earthLuna = CR3BPUtilities.lagrangeStability(mu_earthLuna);
eigenvalue = lagrangeStabilities_earthLuna.InPlaneEigenvalue3(1);

dRDotVec_0 = CR3BPUtilities.librationPlanarOrbitVelocity(rVec_L1, dRVec_0, eigenvalue, mu_earthLuna)
dXVec_0 = [dRVec_0; dRDotVec_0];
 
P = CR3BPUtilities.librationPeriod(eigenvalue);
xVec_0 = [rVec_L1; 0; 0; 0] + [dXVec_0(1:2); 0; dXVec_0(3:4); 0];

[T, X, Phi] = integrateStm(P, xVec_0, mu_earthLuna);

figure
hold on
plot(rVec_L1(1), rVec_L1(2), 'rd', 'DisplayName', 'L_1')
plot(X(:, 1), X(:, 2), 'DisplayName', 'Nonlinear Orbit')
axis equal
xlabel('x')
ylabel('y')
grid on
title('L_1 orbit')
legend

detPhi = zeros(1, length(T));
for i = 1:length(T)
    detPhi(i) = det(Phi(:,:,i));
end
error = abs(1 - detPhi);

figure
plot(T, error)
xlabel('Time')
ylabel('Determinant of \Phi')
title('Determinant of State Transition Matrix over Time')
grid on

function dYdt = cr3bpStm(t, YVec, mu)
    xVec = YVec(1:6);
    xVecDot = CR3BPUtilities.cr3bpEom(t, xVec, mu);
    
    Phi = reshape(YVec(7:end), 6, 6);
    A = CR3BPUtilities.linearizedA(xVec(1:3), mu);
    PhiDot = A*Phi;
    
    dYdt = [xVecDot; PhiDot(:)];
end

function [T, X, Phi] = integrateStm(T, xVec_0, mu)
    tspan = [0, T];
    YVec_0 = [xVec_0; reshape(eye(6), [], 1)];

    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [T, Y] = ode113(@(t,Y) cr3bpStm(t,Y,mu), tspan, YVec_0, options);

    X = Y(:,1:6);
    
    numSteps = length(T);
    Phi = zeros(6, 6, numSteps);
    for i = 1:numSteps
        Phi(:,:,i) = reshape(Y(i, 7:end), 6, 6);
    end
end