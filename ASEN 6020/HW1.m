clear all; clc; close all;
addpath('..')

%% 2
J_BE = @(r, l) sqrt(2*l/(1 + l)) - 1 + sqrt(2*r/(l*(l + r))) - sqrt(2/(l*(l + 1))) + sqrt(2*l/(r*(l + r))) - 1/sqrt(r);
J_H = @(r) sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));

J_compare = @(r, l) J_H(r) - J_BE(r, l); % J < 0 means Hohmann is better, J > 0 means Bi-Elliptic is better

r = 1:.1:100;
l = 1:1:100;

J = NaN(length(r), length(l));
for i = 1:length(r)
    for j = 1:length(l)
        J(i, j) = J_compare(r(i), l(j));
    end
end

figure
contourf(r, l, J', 'ShowText', 'on');
colormap(jet); 
colorbar;
xlabel('r');
ylabel('l');
title('Planar Transfer Optimality, J_{Hohmann} - J_{bi-elliptical}');

%% 4
dogLegPlaneChange1090 = @(di_1, di_3) sqrt(31/11 - 2*sqrt(20/11)*cos(di_1)) + sqrt(31/11 - 2*sqrt(20/11)*cos(di_3)) + 2*sqrt(1/55)*sin((pi/2 - di_1 - di_3)/2);

dis = linspace(0, pi/72, 100);

J_planeChange = NaN(length(dis));
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
title('Dog-Leg Plane Change Optimality, r = 10, \Deltai = 90\circ');

[minValue, minIndex] = min(dJ_planeChange(:));
[row, col] = ind2sub(size(dJ_planeChange), minIndex);
fprintf('Minimum dJ_planeChange: %f at di_1 = %f degrees, di_3 = %f degrees\n', minValue, rad2deg(dis(row)), rad2deg(dis(col)));

%% 5
dDvdeta = @(r, di, eta) di*sqrt(2*r/(r + 1))*sin(eta*di)/sqrt(2*r/(r + 1) + 1 - 2*sqrt(2*r/(r + 1))*cos(eta*di)) - di*2/(r*(r + 1))*sin((1 - eta)*di)/(2/sqrt(r*(r + 1))*sqrt(1 - cos((1 - eta)*di)));

dis = pi/12:pi/12:5*pi/12;
rs = 1:1:100;
etaStar = NaN(length(dis), length(rs));
options = optimoptions('fsolve', 'Display', 'off');
for i = 1:length(dis)
    di = dis(i);
    for j = 1:length(rs)
        r = rs(j);
        etaStar(i, j) = fsolve(@(eta) dDvdeta(r, di, eta), .5, options);
    end
end

figure
hold on
for i = 1:length(dis)
    plot(rs, etaStar(i, :), 'DisplayName', ['\Deltai = ', num2str(rad2deg(dis(i))), '\circ']);
end
hold off
xlabel('r');
ylabel('\eta^*');
title('Optimal \eta* for Different Plane Change Angles');
legend

%% 6
dv_b = @(e) 2*(sqrt(1 - e) - (1 - e));
dv_c = @(e) 2*(1 + e - sqrt(1 + e));
es = 0:.1:1;
J = NaN(length(es), 2);
for i = 1:length(es)
    e = es(i);
    J(i, 1) = dv_b(e);
    J(i, 2) = dv_c(e);
end

figure
hold on
plot(es, J(:, 1), 'DisplayName', 'dv_b');
plot(es, J(:, 2), 'DisplayName', 'dv_c');
hold off
xlabel('e');
ylabel('dv');
legend

%% 7
J_1 = @(v_inf) sqrt(v_inf^2 + 2) - 1;
J_2p = @(v_inf, r) sqrt(v_inf^2 + 2/r) + 1 - sqrt(2)/r;
J_2a = @(v_inf, r) sqrt(v_inf^2 + 2/r) - 1 + sqrt(2)*(r - 1)/sqrt(r*(r + 1));

g_I = @(v_inf, r) J_1(v_inf) - J_2p(v_inf, r);
g_II = @(v_inf, r) J_1(v_inf) - J_2a(v_inf, r);

hyperbolicDepartureOptimal(g_I);
title('g_I = J_1 - J_{2p}');
xlim([0, 1]);

hyperbolicDepartureOptimal(g_II);
title('g_{II} = J_1 - J_{2a}');
xlim([1, 10]);

function hyperbolicDepartureOptimal(g)
    rs = 0:.01:10;
    vs = 0:.01:2;
    J = NaN(length(rs), length(vs));
    for i = 1:length(rs)
        for j = 1:length(vs)
            J(i, j) = g(vs(j), rs(i));
        end
    end

    figure
    contour(rs, vs, J', 'ShowText', 'on', 'LineColor', 'k');
    xlabel('l/r');
    ylabel('v_{\infty}/v_c');
end