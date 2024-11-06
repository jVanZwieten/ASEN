function a_solution = lambertIterator(c, s,  mu, longTOF, longWay, hyperbolic, TOF_known, a_guess, toleranceFactor)
    if(hyperbolic)
        deltaTOF = @(a) lambertsEquationHyperbolic(c, s, a, mu, longWay) - TOF_known;
    else
        deltaTOF = @(a) lambertsEquation(c, s, a, mu, longTOF, longWay) - TOF_known;
    end

    iteratorOptions = optimoptions('fsolve', 'TolFun', TOF_known*toleranceFactor,'Display', 'off');
    
    for delta = .002:.002:1
        [a_solution, fval, exitFlag] = fsolve(deltaTOF, a_guess*(1 + delta), iteratorOptions);
        if(exitFlag == 1) break; end
    end
    assert(exitFlag == 1 && imag(a_solution) == 0)
end

function TOF = lambertsEquationHyperbolic(c, s, a, mu, longWay)
    longWayFactor = -1 + 2*longWay;

    n = sqrt(mu/a^3);
    alpha = 2*asinh(sqrt(s/(2*abs(a))));
    beta = 2*asinh(sqrt((s - c)/(2*abs(a))));
    
    TOF = (sinh(alpha) - alpha + longWayFactor*(sinh(beta) - beta))/n;
end

function [TOF, alpha, beta, n] = lambertsEquation(c, s, a, mu, longTof, longWay)
    n = sqrt(mu/a^3);

    alpha = 2*asin(sqrt(s/(2*a)));
    if(longTof) alpha = 2*pi - alpha; end

    beta = 2*asin(sqrt((s - c)/(2*a)));
    if(longWay) beta = -beta; end
    
    TOF = lambertsEquationGeometric(alpha, beta, n);
end

function TOF = lambertsEquationGeometric(alpha, beta, n)
    TOF = ((alpha - beta) - (sin(alpha) - sin(beta)))/n;
end