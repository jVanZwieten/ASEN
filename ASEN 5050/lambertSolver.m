function [a, e, nu_1, nu_2, vVec_x1, vVec_x2]  = lambertSolver(rVec_1, rVec_2, TOF, longWay, mu)
    AU = astroUtilities;
    toleranceFactor = 1e-10;

    [r_n, rHat_n] = Utilities.magnitudeDirection([rVec_1 rVec_2]);
    r_1 = r_n(1);
    r_2 = r_n(2);
    rHat_1 = rHat_n(:, 1);
    rHat_2 = rHat_n(:, 2);

    deltaNu = Utilities.normalizeAngle(Utilities.angleBetweenUnitVectors(rHat_1, rHat_2));
    if(longWay) deltaNu = 2*pi - deltaNu; end

    [c, s] = lambertGeometricQuantities(r_1, r_2, deltaNu);

    TOF_parabolic = lambertsEquationParabolic(c, s, mu, longWay);
    hyperbolicTransfer = TOF < TOF_parabolic;

    TOF_min = lambertsEquationMin(c, s, longWay, mu);
    longTOF = TOF > TOF_min;
    a_guess = s/2;

    a = lambertIterator(c, s, mu, longTOF, longWay, hyperbolicTransfer, TOF, a_guess, toleranceFactor);

    [TOF_result, alpha, beta] = lambertsEquation(c, s, a, mu, longTOF, longWay);
    assert(abs(TOF - TOF_result)/TOF < toleranceFactor*10, "TOF: %s; Result: %s", TOF, TOF_result)

    p = semiLatusRectum(a, c, s, r_1, r_2, alpha, beta);
    e = AU.eccentricityFromSemiLatusRectum(p, a);

    nu_1 = AU.trueAnomalyFrompre(p, r_1, e);
    nu_2 = AU.trueAnomalyFrompre(p, r_2, e);

    if(~deltaNuCheck(nu_1, nu_2, deltaNu))
        nu_1 = -nu_1;
    end
    if(~deltaNuCheck(nu_1, nu_2, deltaNu))
        nu_2 = -nu_2;
    end
    if(~deltaNuCheck(nu_1, nu_2, deltaNu))
        nu_1 = -nu_1;
    end
    if(~deltaNuCheck(nu_1, nu_2, deltaNu))
        error("Couldn't find a combination of true anomalies which satisfied the expected deltaNu: \n deltaNu: %s \n nu1: %s \n nu2: %s", deltaNu, abs(nu_1), abs(nu_2))
    end

    f = AU.fFunction(r_2, p, deltaNu);
    g = AU.gFunction(r_2, r_1, mu, p, deltaNu);
    vVec_x1 = (rVec_2 - f*rVec_1)/g;

    fDot = AU.fDotFunction(mu, p, deltaNu, r_2, r_1);
    gDot = AU.gDotFunction(r_1, p, deltaNu);
    vVec_x2 = fDot*rVec_1 + gDot*vVec_x1;
end

function [c, s] = lambertGeometricQuantities(r_1, r_2, deltaNu)
    c = sqrt(r_1^2 + r_2^2 - 2*r_1*r_2*cos(deltaNu));
    s = (r_1 + r_2 + c)/2;
end

function TOF_parabolic = lambertsEquationParabolic(c, s, mu, longWay)
    longWayFactor = -1 + 2*longWay;
    TOF_parabolic = (sqrt(2/mu)*(s^(3/2) + longWayFactor*(s - c)^(3/2)))/3;
end

function TOF_min = lambertsEquationMin(c, s, longWay, mu)
    a_min = s/2;
    n_min = sqrt(mu/a_min^3);
    alpha_min = pi;
    beta_min = 2*asin(sqrt((s - c)/s));
    if(longWay) beta_min = -beta_min; end
    
    TOF_min = lambertsEquationGeometric(alpha_min, beta_min, n_min);
end

function TOF = lambertsEquationGeometric(alpha, beta, n)
    TOF = ((alpha - beta) - (sin(alpha) - sin(beta)))/n;
end

function [TOF, alpha, beta, n] = lambertsEquation(c, s, a, mu, longTof, longWay)
    n = sqrt(mu/a^3);

    alpha = 2*asin(sqrt(s/(2*a)));
    if(longTof) alpha = 2*pi - alpha; end

    beta = 2*asin(sqrt((s - c)/(2*a)));
    if(longWay) beta = -beta; end
    
    TOF = lambertUtilities.lambertsEquationGeometric(alpha, beta, n);
end

function p = semiLatusRectum(a, c, s, r_1, r_2, alpha, beta)
    p = 4*a*(s - r_1)*(s - r_2)/c^2*sin((alpha + beta)/2)^2;
end

function result = deltaNuCheck(nu_1, nu_2, deltaNu)
    result = abs(nu_2 - nu_1 - deltaNu) < .1;
end