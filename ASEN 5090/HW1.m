addpath("..")
clear

% 1-2
% a
v_0=50; % m/s
x_0 = 250; % m
h_0 = 100; % m

x = -x_0:v_0:x_0;
T = 0:length(x)-1;
X = [x;
    h_0 * ones(1, length(x));];

Y = [abs(atan2(X(1, :), X(2, :)));
    sqrt(X(1, :).^2 + X(2, :).^2)];
Y = [Y;
    v_0*sin(pi - Y(1, :))];

YDeg = [rad2deg(Y(1, :));
    Y(2:3, :)];

titles = ["zenith angle (degrees)" "range (m)" "range rate (m/s)"];
figure
tiledlayout(3, 1)
for i = 1:size(Y, 1)
    nexttile
    plot(T, YDeg(i, :))
    title(titles(i))
end

% 2-1
% a
bits = 1023;
[CA_19, G_19_1, G_19_2] = GnssUtilities.generateCA([3 6], bits);
plotCAFirstLast16(CA_19, '19')

% b
CA_19_b = GnssUtilities.generateCAFromRegisters([3 6], G_19_1, G_19_2, bits);
plotCAFirstLast16(CA_19_b, '19 (bits 1024-2046)')

% c
CA_25 = GnssUtilities.generateCA([5 7], bits);
plotCAFirstLast16(CA_25, '25')

% d
CA_5 = GnssUtilities.generateCA([1 9], bits);
plotCAFirstLast16(CA_5, '5')

% 2-2
% a
R19N = GnssUtilities.autoCorrelateFunction(CA_19);
plotCorrelation(R19N)
title("R^{(19)}(n)")

% b
CA_19_delay200 = Utilities.rightShift(CA_19, 200);
R19Delay200N = GnssUtilities.crossCorrelateFunction(CA_19, CA_19_delay200);
plotCorrelation(R19Delay200N)
title("R^{(19, 19 delayed 200)}(n)")
% The peak is found at 200 chips, which is expected for a signal delayed by 200 chips.

% c
R19_25 = GnssUtilities.crossCorrelateFunction(CA_19, CA_25);
plotCorrelation(R19_25)
title("R^{(19, 25)}(n)")
% This is very similar to the plot in part a, but without the spike to R = 1. This is of course because the signals are not just out of phase, they are different signals.

% d
R19_5 = GnssUtilities.crossCorrelateFunction(CA_19, CA_5);
plotCorrelation(R19_5)
title("R^{(19, 5)}(n)")
% This is essentially the same as part c, for the same reason: the signals are fundamentally different, not just unsynced.

% e
x1 = GnssUtilities.CaToX(Utilities.rightShift(CA_19, 350));
x2 = GnssUtilities.CaToX(Utilities.rightShift(CA_25, 905));
x3 = GnssUtilities.CaToX(Utilities.rightShift(CA_5, 75));
xSum = x1 + x2 + x3;
x19 = GnssUtilities.CaToX(CA_19);
R19_xSum = GnssUtilities.crossCorrelateFunctionX(x19, xSum);
plotCorrelation(R19_xSum)
title("R^{(19, 19 + 25 + 5)}(n)")
% The peak is at n=350. This makes sense, given at n=350, CA_19 is synced with the delayed CA_19 signal, even though there's no correlation to the other two signals mixed in.

% f
noise = 4 * randn(1, 1023);
signals = {x1, x2, x3, noise};
T = 1:1023;
titles=["x1" "x2" "x3" "noise"];
figure
tiledlayout(4, 1)
for i = 1:length(signals)
    nexttile
    plot(T, signals{i})
    ylim([-15 15])
    xlim([0 1023])
    title(titles(i))
end

% g
noisySum = x1 + x2 + x3 + noise;
R19_noisySum = GnssUtilities.crossCorrelateFunctionX(x19, noisySum);
plotCorrelation(R19_noisySum)
title("R^{(19, 19 + 25 + 5 + noise)}(n)")
%% The correlation spikes at n = 350, correctly identifying the start of PRN19's CA. I honestly am surprised that this algorithm so clearly identifies the start of the CA code that it's seeking, but it does explain why GPS is so effective despite the broadcast being so weak, even below the noise floor.

function plotCAFirstLast16(G, prn)
    first16 = G(1:16);
    firstHex = Utilities.binVector2Hex(first16);
    last16 = G(1008:1023);
    lastHex = Utilities.binVector2Hex(last16);

    figure
    tiledlayout(2, 1)
    nexttile
    plot(1:16, first16, 'o')
    title(['PRN ' prn ' chips 1-16 (' firstHex ')'])
    nexttile
    plot(1008:1023, last16, 'o')
    title(['PRN ' prn ' chips 1008-1023 (' lastHex ')'])
end

function plotCorrelation(R)
    figure
    plot(R(1, :), R(2, :))
end