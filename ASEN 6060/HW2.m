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

%% 1c
muFactors = [.5, 2];
for muFactor = muFactors
    mu = mu_earthLuna*muFactor;
    lagrangePoints = CR3BPUtilities.lagrangePoints(mu);
    CR3BPUtilities.lagrangePlot(lagrangePoints)
    title(['Lagrange Points, \mu factor =', num2str(muFactor)])
end

%% 2
lagrangeStabilities_earthLuna = CR3BPUtilities.langrangeStability(mu_earthLuna);
lagrangeStabilities_solEarth = CR3BPUtilities.langrangeStability(CR3BPUtilities.mu_solEarth);