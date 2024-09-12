% % 1-2
% v_0=50; % m/s
% x_0 = 250; % m
% h_0 = 100; % m

% x = -x_0:v_0:x_0;
% T = 0:1:length(x)-1;
% X = [x;
%     h_0 * ones(1, length(x));];

% Y = [atan2(X(1, :), X(2, :));
%     sqrt(X(1, :).^2 + X(2, :).^2)];
% Y = [Y(1:2, :);
%     v_0*sin(pi - Y(1, :))];

% YDeg = [rad2deg(Y(1, :));
%     Y(2:3, :)];

% figure
% plot(T', YDeg', 'o')
% legend({"zenith angle (degrees)" "range (m)" "range rate (m/s)"})

% 2-1
% a
iterations = 1023;
G1 = ones(1, 10);
G2 = ones(1, 10);
G = NaN(1, iterations);

for i = 1:1:iterations
    G1i = G1(10);
    
    S = [G2(3) G2(6)];
    G2i = mod(sum(S), 2);

    XGi = mod(sum([G1i G2i]), 2);
    G(i) = XGi;

    G1Feedback = [G1(3) G1(10)];
    G1 = [mod(sum(G1Feedback), 2) G1(1:9)];

    G2Feedback = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2 = [mod(sum(G2Feedback), 2) G2(1:9)];
end

G