function X = ballisticTrajectorySolver(dt, t_f, v_0, gamma_0, h_0, m, dragCoefficient, area, simplified)
    HU = HypersonicsUtilities;
    g = HU.earthSurfaceGravityAcceleration;
    w = m*g;

    steps = t_f/dt + 1; % +1 accounts for t = 0;
    X = [zeros(5, length(steps))]; % [t, v, gamma, r, theta]
    X(:, 1) = [0; v_0; gamma_0; h_0 + HU.earthRadius; 0];
    
    for i = 2:steps
        X_prior = X(:, i - 1);

        h = HU.altitude(X(4, i - 1));
        rho = HU.airDensity(h);
        drag = HU.drag(dragCoefficient, area, rho, X_prior(2));

        if simplified
            vDot = HU.decelerationBallisticSimple(m, drag);
            gammaDot = 0;
        else
            vDot = HU.acceleration(0, drag, m, X_prior(2));
            gammaDot = HU.flightPathAngleRate(X_prior(2), 0, w, X_prior(4), X(3));
        end

        A = [1
        vDot
        gammaDot
        HU.radialDistanceRate(X_prior(2), X_prior(3))
        HU.earthAngleRate(X_prior(2), X_prior(3), X_prior(4))];

        X(:, i) = EulerTimeStep(X_prior, dt, A);
    end
end

function X_next = EulerTimeStep(X, dt, A)
    X_next = X + dt*A;
end