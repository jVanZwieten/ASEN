function dYdt = cr3bpStm(t, YVec, mu)
    xVec = YVec(1:6);
    xVecDot = CR3BPUtilities.cr3bpEom(t, xVec, mu);
    
    Phi = reshape(YVec(7:end), 6, 6);
    A = compute_A_matrix(xVec, mu);
    PhiDot = A*Phi;
    
    dYdt = [xVecDot; PhiDot(:)];
end

function Phi_tf_t0 = integrateCr3bp(T, xVec_0, mu)
    tspan = [0, T];
    YVec_0 = [xVec_0; reshape(eye(6), [], 1)];

    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [t, Y] = ode113(@(t,Y) cr3bpStm(t,Y,mu), tspan, Y0, options);

    Phi_tf_t0 = reshape(Y(end,7:end), 6, 6);
end