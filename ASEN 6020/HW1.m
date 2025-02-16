clear all; clc; close all;

% Define the cost functions for Bi-Elliptic and Hohmann transfers
J_BE = @(r, l) sqrt(2*l/(1 + l)) - 1 + sqrt(2*r/(l*(l + r))) - sqrt(2/(l*(l + 1))) + sqrt(2*l/(r*(l + r))) - 1/sqrt(r);
J_H = @(r) sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));

% Define the comparison function
J_compare = @(r, l) J_H(r) - J_BE(r, l); % J < 0 means Hohmann is better, J > 0 means Bi-Elliptic is better

% Define the range for r and l
r = 1:.1:20; % Range for r
l = 1:1:20; % Range for l

% Initialize the matrix to store the comparison results
J = zeros(length(r), length(l));

% Compute the comparison results
for i = 1:length(r)
    for j = 1:length(l)
        J(i, j) = J_compare(r(i), l(j));
    end
end

figure
% Plot the results using a contour plot
contourf(r, l, J', 'ShowText', 'on'); % Note the transpose of J for correct orientation
colormap(jet); % You can choose other colormaps like 'parula', 'hsv', etc.
colorbar; % Adds a colorbar to the plot
xlabel('r');
ylabel('l');
title('Comparison of Hohmann and Bi-Elliptic Transfers');

%% 4
dogLegPlaneChange1090 = @(di_1, di_3) sqrt(31/11 - 2*sqrt(20/11)*cos(di_1)) + sqrt(31/11 - 2*sqrt(20/11)*cos(di_3)) + 2*sqrt(1/55)*sin((pi/2 - di_1 - di_3)/2);

dis = linspace(0, pi/72, 100);

J_planeChange = NaN(length(dis), length(dis));
for i = 1:length(dis)
    for j = 1:length(dis)
        if(dis(i) + dis(j) < pi/2)
            J_planeChange(i, j) = dogLegPlaneChange1090(dis(i), dis(j));
        end
    end
end

dv_ib = 2*(sqrt(2*10/11) - 1 + sqrt(2/110)*sin(pi/4));
Utilities.assertIsWithin(dv_ib, J_planeChange(1, 1), 1e-6);
dJ_planeChange = J_planeChange - dv_ib;

figure
contourf(rad2deg(dis), rad2deg(dis), dJ_planeChange', 'ShowText', 'on');
colormap(jet);
colorbar;
xlabel('\Deltai_1 (\circ)');
ylabel('\Deltai_3 (\circ)');
title('Dog-Leg Plane Change, r = 10, \Deltai = 90\circ');

% Find the minimum value in dJ_planeChange and its indices
[minValue, minIndex] = min(dJ_planeChange(:));
[row, col] = ind2sub(size(dJ_planeChange), minIndex);

% Display the minimum value and corresponding angles
fprintf('Minimum dJ_planeChange: %f at di_1 = %f degrees, di_3 = %f degrees\n', minValue, rad2deg(dis(row)), rad2deg(dis(col)));