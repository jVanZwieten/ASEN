addpath('..')
addpath('..\..')

clear
close all

alphasA = [0 pi/6 pi/3];

flatPlateFileName = "Flat_Plate_Processed.dat";
flatPlateReferenceArea = 1; % m^2
[body_flatPlate, aeroCoefficients_FlatPlate] = NewtonianAerodynamics(flatPlateFileName, alphasA, flatPlateReferenceArea);
[theoreticalC_L, theoreticalC_D] = plateAeroforceCoefficients(alphasA);
theoreticalC_L = theoreticalC_L';
theoreticalC_D = theoreticalC_D';
alpha = alphasA';
theoreticalAeroCoefficients_FlatPlate = table(rad2deg(alpha), theoreticalC_L, theoreticalC_D);

sphereFileName = "Sphere_Processed.dat";
sphereReferenceArea = pi; % m^2
[body_sphere, aeroCoefficients_sphere] = NewtonianAerodynamics(sphereFileName, alphasA, sphereReferenceArea);

apolloFileName = "Apollo_Processed.dat";
apolloReferenceArea = 1.9558^2*pi; % m^2
alphasA = deg2rad(-20);
[body_apollo, aeroCoefficients_apollo] = NewtonianAerodynamics(apolloFileName, alphasA, apolloReferenceArea);
LDRatio_apollo = aeroCoefficients_apollo.C_L(1)/aeroCoefficients_apollo.C_D(1);

function [C_L, C_D] = plateAeroforceCoefficients(alpha)
    C_D = 2*sin(alpha).^3;
    C_L = 2*sin(alpha).^2.*cos(alpha);
end