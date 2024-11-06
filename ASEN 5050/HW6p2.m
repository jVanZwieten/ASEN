clear
close all
addpath('..')

mu_sun = CelestialParameters.gravityParameter_sun;

load("HW6EphemEarth.mat")
load("HW6EphemMars.mat")

varNames = {'Epoch_Earth', 'Epoch_Mars', 'v_infinityEarth', 'v_infinityMars'};
varTypes = {'double', 'double', 'double', 'double'};
transfersT = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
for i_earth = 1:size(HW6EphemEarth, 1)
    x_earthI = HW6EphemEarth(i_earth, :);
    rVec_1 = [x_earthI.X; x_earthI.Y; x_earthI.Z];
    t_1 = x_earthI.Epoch;

    x_mars = HW6EphemMars(HW6EphemMars.Epoch > t_1, :);
    for i_mars = 1:size(x_mars, 1)
        x_marsI = x_mars(i_mars, :);
        rVec_2 = [x_marsI.X; x_marsI.Y; x_marsI.Z];
        t_2 = x_marsI.Epoch;

        TOF = Utilities.days2s(t_2 - t_1);
        [~, ~, ~, ~, vVec_x1, vVec_x2] = lambertSolver(rVec_1, rVec_2, TOF, false, mu_sun);
        transfersT = [transfersT; {t_1, t_2, norm(vVec_x1), norm(vVec_x2)}];
    end
end