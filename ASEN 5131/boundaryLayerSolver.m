function [eta, X] = boundaryLayerSolver(M_e, b_0, g_0, e_0)
    Pr = 0.75;
    gamma = 1.4;
    n = 0.72;

    eta_range = [0 20];
    X_0 = [0; 0; b_0; g_0; e_0; 0]; % f, a, b, g, e, z

    ode_function = @(eta, X) boundaryLayerODE(eta, X, Pr, gamma, M_e, n);
    [eta, X] = ode45(ode_function, eta_range, X_0);

    function dydeta = boundaryLayerODE(eta, X, Pr, gamma, M_e, n)
        % X = [f, a, b, g, e, z]
        f = X(1);  % f
        a = X(2);  % f'
        b = X(3);  % f''
        g = X(4);  % g = h/h_e
        e = X(5);  % g'
        z = X(6);  % y√(Re_x)/x
        
        C = g^(n - 1);
        Cprime = (n - 1) * g^(n - 2) * e;
        
        b_prime = -(b / C) * (Cprime + f); % Eq6
        e_prime = -(Cprime / C) * e - Pr / C * f * e - Pr * (gamma - 1) * M_e^2 * b^2; % Eq7
        z_prime = sqrt(2)*g; % New ODE for y√(Re_x)/x
        
        dydeta = [a; b; b_prime; e; e_prime; z_prime];
    end
end