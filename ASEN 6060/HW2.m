clc; clc; clear; close all;
format long e
addpath('..')

%% 1a
mu_earthLuna = CR3BPUtilities.mu_earthLuna
lagrangePoints_earthLuna = CR3BPUtilities.lagrangePoints(mu_earthLuna);

for i = 1:5
    lagrangePoints_earthLuna{i}
end

CR3BPUtilities.lagrangePlot(lagrangePoints_earthLuna)
title("Earth-Lunar Lagrange Points")

%% 1b
CM = NaN(1, 5);
for i = 1:length(CM)
    CM(i) = CR3BPUtilities.jacobiConstant(lagrangePoints_earthLuna{i}, mu_earthLuna);
end

CM

%% 1c
figure
muFactors = [.5, 2];
for muFactor = muFactors
    mu = mu_earthLuna*muFactor;
    lagrangePoints = CR3BPUtilities.lagrangePoints(mu);
    CR3BPUtilities.lagrangePlot(lagrangePoints)
    title(['Lagrange Points, \mu factor =', num2str(muFactor)])
end

%% 2
lagrangeStabilities_earthLuna = CR3BPUtilities.lagrangeStability(mu_earthLuna);
lagrangeStabilities_solEarth = CR3BPUtilities.lagrangeStability(CR3BPUtilities.mu_solEarth);

%% 4c
rVec_L1 = lagrangePoints_earthLuna{1};
dRVec_0 = [-.001; 0];
eigenvalue = lagrangeStabilities_earthLuna.InPlaneEigenvalue3(1);

dRDotVec_0 = CR3BPUtilities.librationPlanarOrbitVelocity(rVec_L1, dRVec_0, eigenvalue, mu_earthLuna)
dXVec_0 = [dRVec_0; dRDotVec_0];

P = CR3BPUtilities.librationPeriod(eigenvalue);
tSpan = [0 P]

[t, X_linear] = CR3BPUtilities.integrateLinearizedCr3bp(rVec_L1, dXVec_0, mu_earthLuna, tSpan, -14);
X_linear(:, 1) = X_linear(:, 1) + rVec_L1(1);
X_linear(:, 2) = X_linear(:, 2) + rVec_L1(2);

figure
hold on
plot(rVec_L1(1), rVec_L1(2), 'rd', 'DisplayName', 'L_1')
plot(X_linear(:, 1), X_linear(:, 2), 'DisplayName', 'Linearized Orbit')
axis equal
xlabel('x')
ylabel('y')
grid on
title('L_1 orbit')

xVec_0 = [rVec_L1; 0; 0; 0] + [dXVec_0(1:2); 0; dXVec_0(3:4); 0];
[t, X_nonLinear] = CR3BPUtilities.integrateCr3bp(xVec_0, mu_earthLuna, tSpan, -14);
plot(X_nonLinear(:, 1), X_nonLinear(:, 2), 'DisplayName', 'Nonlinear Orbit')
legend

%% 4e
A2D_earthLunaL1 = CR3BPUtilities.equilibriumPerturbationA2D(rVec_L1, mu_earthLuna);
[V, D] = eig(A2D_earthLunaL1);
lambda_3 = D(3, 3)
lambda_4 = D(4, 4)
eigenvector_3 = V(:, 3)
eigenvector_4 = V(:, 4)

basisVec_1 = real(eigenvector_3)
basisVec_2 = imag(eigenvector_3)
dXhat_0 = Utilities.unitVector(dXVec_0)
basisHat_1 = Utilities.unitVector(basisVec_1)
dXVec_0(1)/basisVec_1(1)*basisVec_1

%% 5a
rVec_L4 = lagrangePoints_earthLuna{4};
A2D_earthLunaL4 = CR3BPUtilities.equilibriumPerturbationA2D(rVec_L4, mu_earthLuna);
[eigenvectors_L4, eigenvalues_L4] = eig(A2D_earthLunaL4)
eigenvalues_L4 = j*imag(eigenvalues_L4*ones(4, 1))

P_1 = CR3BPUtilities.librationPeriod(eigenvalues_L4(1))
P_3 = CR3BPUtilities.librationPeriod(eigenvalues_L4(3))

eigenvector_L4short = eigenvectors_L4(:, 1);
basisVec_1L4short = real(eigenvector_L4short);
dXVec_0L4short = .02/norm(basisVec_1L4short(1:2))*basisVec_1L4short
norm(dXVec_0L4short(1:2))

[T, X_linearL4short] = CR3BPUtilities.integrateLinearizedCr3bp(rVec_L4, dXVec_0L4short, mu_earthLuna, [0 P_1], -14);
X_linearL4short(:, 1) = X_linearL4short(:, 1) + rVec_L4(1);
X_linearL4short(:, 2) = X_linearL4short(:, 2) + rVec_L4(2);
xVec_0L4short = [rVec_L4; 0; 0; 0] + [dXVec_0L4short(1:2); 0; dXVec_0L4short(3:4); 0];
[T, X_nonLinearL4short] = CR3BPUtilities.integrateCr3bp(xVec_0L4short, mu_earthLuna, [0 P_1], -14);
figure
hold on
plot(rVec_L4(1), rVec_L4(2), 'rd', 'DisplayName', 'L_1')
plot(X_linearL4short(:, 1), X_linearL4short(:, 2), 'DisplayName', 'Linearized Orbit')
plot(X_nonLinearL4short(:, 1), X_nonLinearL4short(:, 2), 'DisplayName', 'Nonlinear Orbit')
plot(xVec_0L4short(1), xVec_0L4short(2), 'ro', 'DisplayName', '$\vec{r}_0$')
plot(X_nonLinearL4short(end, 1), X_nonLinearL4short(end, 2), 'o', 'DisplayName', 'Nonlinearized $\vec{r}_f$')
axis equal
xlabel('x')
ylabel('y')
legend('Interpreter', 'latex')
grid on
title('L4 short orbit')

eigenvector_L4long = eigenvectors_L4(:, 3);
basisVec_1L4long = real(eigenvector_L4long);
dXVec_0L4long = .02/norm(basisVec_1L4long(1:2))*basisVec_1L4long
norm(dXVec_0L4long(1:2))

[T, X_linearL4long] = CR3BPUtilities.integrateLinearizedCr3bp(rVec_L4, dXVec_0L4long, mu_earthLuna, [0 P_3], -14);
X_linearL4long(:, 1) = X_linearL4long(:, 1) + rVec_L4(1);
X_linearL4long(:, 2) = X_linearL4long(:, 2) + rVec_L4(2);
xVec_0L4long = [rVec_L4; 0; 0; 0] + [dXVec_0L4long(1:2); 0; dXVec_0L4long(3:4); 0];
[T, X_nonLinearL4long] = CR3BPUtilities.integrateCr3bp(xVec_0L4long, mu_earthLuna, [0 P_3], -14);
figure
hold on
plot(rVec_L4(1), rVec_L4(2), 'rd', 'DisplayName', 'L_1')
plot(X_linearL4long(:, 1), X_linearL4long(:, 2), 'DisplayName', 'Linearized Orbit')
plot(X_nonLinearL4long(:, 1), X_nonLinearL4long(:, 2), 'DisplayName', 'Nonlinear Orbit')
plot(xVec_0L4long(1), xVec_0L4long(2), 'ro', 'DisplayName', '$\vec{r}_0$')
plot(X_nonLinearL4long(end, 1), X_nonLinearL4long(end, 2), 'o', 'DisplayName', 'Nonlinearized $\vec{r}_f$')
axis equal
xlabel('x')
ylabel('y')
legend('Interpreter', 'latex')
grid on
title('L_4 long orbit')